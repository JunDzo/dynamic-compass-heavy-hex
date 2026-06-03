#!/usr/bin/env python3
"""
solve_surface.py

Library for solving and plotting fitted surface parameters.

Supports:
- Linear: z = a + b x + c y
- Quadratic: z = a + b x + c y + d x^2 + e x y + f y^2

Can be used as a library or as a CLI tool.

CLI Usage examples:
  # Solve z given x and y using linear model if present, else quadratic
  python solve_surface.py --model fit_X_b1p_b2p_p_th.json --unknown z --x 0.02 --y 0.1

  # Solve y given x and z
  python solve_surface.py --model fit_X_b1p_b2p_p_th.json --unknown y --x 0.02 --z 3e-3

  # Force quadratic even if linear exists
  python solve_surface.py --model fit_X_b1p_b2p_p_th.json --unknown x --y 0.07 --z 4e-3 --use quadratic
"""

import argparse
import json
import math
import sys
from typing import List, Tuple, Optional
import matplotlib.pyplot as plt  # plotting

def clamp01(u: float) -> bool:
    """Check if value is in [0, 1] range."""
    return 0.0 <= u <= 1.0


def _solve_linear(unknown: str, a: float, b: float, c: float, x: float=None, y: float=None, z: float=None) -> List[float]:
    """Solve linear equation z = a + b*x + c*y for the unknown variable."""
    # z = a + b x + c y
    if unknown == "z":
        return [a + b*x + c*y]
    if unknown == "y":
        if abs(c) < 1e-300:
            raise ZeroDivisionError("Linear solve for y failed: c ≈ 0.")
        return [(z - a - b*x)/c]
    if unknown == "x":
        if abs(b) < 1e-300:
            raise ZeroDivisionError("Linear solve for x failed: b ≈ 0.")
        return [(z - a - c*y)/b]
    raise ValueError("unknown must be one of {x,y,z}")

def _real_quadratic_roots(A: float, B: float, C: float) -> List[float]:
    """Solve quadratic equation A*u^2 + B*u + C = 0."""
    # A u^2 + B u + C = 0
    if abs(A) < 1e-300:
        if abs(B) < 1e-300:
            return []  # Degenerate
        return [-C / B]
    D = B*B - 4*A*C
    if D < 0:
        return []
    sqrtD = math.sqrt(D)
    return [(-B - sqrtD)/(2*A), (-B + sqrtD)/(2*A)]


def _prob_quadratic_roots(A: float, B: float, C: float) -> List[float]:
    """Solve quadratic and return only roots in [0, 1]."""
    roots = _real_quadratic_roots(A, B, C)
    return [r for r in roots if clamp01(r)]


def _solve_quadratic(unknown: str, coef: List[float], x: float=None, y: float=None, z: float=None) -> List[float]:
    """Solve quadratic surface equation z = a + b*x + c*y + d*x^2 + e*x*y + f*y^2 for unknown."""
    a,b,c,d,e,f = coef
    if unknown == "z":
        zval = a + b*x + c*y + d*x*x + e*x*y + f*y*y
        return [zval] if clamp01(zval) else []
    if unknown == "y":
        A = f
        B = c + e*x
        C = a + b*x + d*x*x - z
        return _prob_quadratic_roots(A, B, C)
    if unknown == "x":
        A = d
        B = b + e*y
        C = a + c*y + f*y*y - z
        return _prob_quadratic_roots(A, B, C)
    raise ValueError("unknown must be one of {x,y,z}")

def find_curve(axis: str, coef: List[float], x: float=None, y: float=None, z: float=None, num: int = 200, monotone: str = "decreasing") -> List[Tuple[float, float]]:
    """
    Generate curve points for a 2D slice through the fitted surface.

    Args:
        axis: One of {'xy', 'xz', 'yz', 'yx', 'zx', 'zy'}
        coef: Coefficient list (length 3 for linear, 6 for quadratic)
        x, y, z: Fixed coordinate values (provide appropriate ones based on axis)
        num: Number of sample points along the curve
        monotone: One of {'none', 'decreasing', 'increasing'} - restricts to monotone intervals

    Returns:
        List of (x, y) tuples representing the curve points
    """
    if axis not in {"xy", "xz", "yz", "yx", "zx", "zy"}:
        raise ValueError("axis must be one of {'xy','xz','yz','yx','zx','zy'}")
    axis = ''.join(sorted(axis))

    if len(coef) == 3:
        model = "linear"
        a, b, c = map(float, coef)
    elif len(coef) == 6:
        model = "quadratic"
        a, b, c, d, e, f = map(float, coef)
    else:
        raise ValueError("coef must have length 3 (linear) or 6 (quadratic)") 

    pts: List[Tuple[float, float]] = []

    # Helper to restrict domain for monotonic behavior
    def restrict_interval_for_monotone(axis_norm: str) -> Tuple[float, float]:
        # Returns [lo, hi] for the sampled variable
        if monotone == "none":
            return 0.0, 1.0
        if model == "linear":
            if axis_norm == 'yz':
                slope = c
            elif axis_norm == 'xz':
                slope = b
            else:
                return 0.0, 1.0
            if monotone == "decreasing" and slope <= 0:
                return 0.0, 1.0
            if monotone == "increasing" and slope >= 0:
                return 0.0, 1.0
            # otherwise no interval satisfies requested monotonicity
            return 0.0, -1.0  # empty
        # quadratic
        if axis_norm == 'yz':
            if x is None:
                return 0.0, -1.0
            # dz/dy = c + e x + 2 f y
            A1 = 2*f
            B1 = c + e*x
            if abs(A1) < 1e-300:
                slope = B1
                if (monotone == "decreasing" and slope <= 0) or (monotone == "increasing" and slope >= 0):
                    return 0.0, 1.0
                return 0.0, -1.0
            y_thresh = -B1/(A1)
            if monotone == "decreasing":
                if A1 > 0:
                    lo, hi = 0.0, min(1.0, y_thresh)
                else:
                    lo, hi = max(0.0, y_thresh), 1.0
            else:  # increasing
                if A1 > 0:
                    lo, hi = max(0.0, y_thresh), 1.0
                else:
                    lo, hi = 0.0, min(1.0, y_thresh)
            return lo, hi
        if axis_norm == 'xz':
            if y is None:
                return 0.0, -1.0
            # dz/dx = b + 2 d x + e y
            A1 = 2*d
            B1 = b + e*y
            if abs(A1) < 1e-300:
                slope = B1
                if (monotone == "decreasing" and slope <= 0) or (monotone == "increasing" and slope >= 0):
                    return 0.0, 1.0
                return 0.0, -1.0
            x_thresh = -B1/(A1)
            if monotone == "decreasing":
                if A1 > 0:
                    lo, hi = 0.0, min(1.0, x_thresh)
                else:
                    lo, hi = max(0.0, x_thresh), 1.0
            else:  # increasing
                if A1 > 0:
                    lo, hi = max(0.0, x_thresh), 1.0
                else:
                    lo, hi = 0.0, min(1.0, x_thresh)
            return lo, hi
        return 0.0, 1.0

    if axis == "xy":
        if z is None:
            raise ValueError("For axis='xy' you must provide z")
        lo, hi = 0.0, 1.0
        for i in range(num):
            xi = lo + (hi - lo) * i/(num-1)
            if model == "linear":
                if abs(c) < 1e-300:
                    continue
                yi = (z - a - b*xi)/c
                if clamp01(yi):
                    pts.append((xi, yi))
            else:
                roots = _solve_quadratic("y", [a,b,c,d,e,f], x=xi, z=z)
                for yi in roots:
                    if clamp01(yi):
                        pts.append((xi, yi))
        return pts

    if axis == "xz":
        if y is None:
            raise ValueError("For axis='xz' you must provide y")
        lo, hi = restrict_interval_for_monotone('xz')
        if hi < lo:
            return []
        for i in range(num):
            xi = lo + (hi - lo) * i/(num-1)
            zi = a + b*xi + c*y if model == "linear" else a + b*xi + c*y + d*xi*xi + e*xi*y + f*y*y
            if clamp01(zi):
                pts.append((xi, zi))
        return pts

    # axis == 'yz'
    if x is None:
        raise ValueError("For axis='yz' you must provide x")
    lo, hi = restrict_interval_for_monotone('yz')
    if hi < lo:
        return []
    for i in range(num):
        yi = lo + (hi - lo) * i/(num-1)
        zi = a + b*x + c*yi if model == "linear" else a + b*x + c*yi + d*x*x + e*x*yi + f*yi*yi
        if clamp01(zi):
            pts.append((yi, zi))
    return pts

def load_model_from_json(json_path: str, use: str = "auto") -> Tuple[str, List[float]]:
    """
    Load fitted surface model from JSON file.

    Args:
        json_path: Path to JSON file with fitted parameters
        use: One of {'auto', 'linear', 'quadratic'}

    Returns:
        Tuple of (model_type, coefficients)

    Raises:
        ValueError if model not found or invalid
    """
    with open(json_path, "r", encoding="utf-8") as f:
        M = json.load(f)

    have_lin = "linear" in M and M["linear"] is not None
    have_quad = "quadratic" in M and M["quadratic"] is not None

    # Choose model
    if use == "linear":
        if not have_lin:
            raise ValueError("Requested linear model, but JSON has no 'linear' section.")
        return ("linear", M["linear"]["coef"])
    elif use == "quadratic":
        if not have_quad:
            raise ValueError("Requested quadratic model, but JSON has no 'quadratic' section.")
        return ("quadratic", M["quadratic"]["coef"])
    else:
        # auto preference: linear first, else quadratic
        if have_lin:
            return ("linear", M["linear"]["coef"])
        elif have_quad:
            return ("quadratic", M["quadratic"]["coef"])
        else:
            raise ValueError("No usable model found in JSON.")


def main():
    p = argparse.ArgumentParser(description="Solve for unknown using fitted surface parameters")
    p.add_argument("--model", required=True, help="Path to JSON file produced by --fit-out")
    p.add_argument("--unknown", choices=["x","y","z"], required=False)
    p.add_argument("--x", type=float, default=None)
    p.add_argument("--y", type=float, default=None)
    p.add_argument("--z", type=float, default=None)
    p.add_argument("--use", choices=["auto","linear","quadratic"], default="auto",
                   help="Which model to use if both are present")

    plot = p.add_argument_group("Plotting")
    plot.add_argument("--plot-curve", action="store_true", help="Plot a 2D curve slice of the surface")
    plot.add_argument("--axis", choices=["xy","xz","yz","yx","zx","zy"], help="2D plane to plot in; provide the fixed third coordinate via --x/--y/--z")
    plot.add_argument("--num", type=int, default=200, help="Number of samples along the free axis")
    plot.add_argument("--out", default=None, help="Save plot as PNG to this path instead of showing interactively")
    plot.add_argument("--monotone", choices=["none","decreasing","increasing"], default="none", help="Restrict the slice to intervals where the curve is monotone in the plotted axis")

    args = p.parse_args()

    try:
        model = load_model_from_json(args.model, args.use)
    except ValueError as e:
        print(f"Error loading model: {e}", file=sys.stderr)
        sys.exit(2)

    # Plotting mode
    if args.plot_curve:
        if not args.axis:
            print("--plot-curve requires --axis {xy,xz,yz} and one fixed coordinate via --x/--y/--z", file=sys.stderr)
            sys.exit(2)
        coef = list(map(float, model[1]))

        try:
            pts = find_curve(args.axis, coef, x=args.x, y=args.y, z=args.z, num=args.num, monotone=args.monotone)
        except ValueError as e:
            print(f"Plot failure: {e}", file=sys.stderr)
            sys.exit(2)

        if not pts:
            print("No points within [0,1] to plot for the specified slice.", file=sys.stderr)
            sys.exit(1)

        xs, ys = zip(*pts)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xs, ys)
        # normalize label order
        axis_norm = ''.join(sorted(args.axis))
        if axis_norm == 'xy':
            ax.set_xlabel('x')
            ax.set_ylabel('y')
        elif axis_norm == 'xz':
            ax.set_xlabel('x')
            ax.set_ylabel('z')
        else:  # 'yz'
            ax.set_xlabel('y')
            ax.set_ylabel('z')
        # Set axis limits dynamically based on data with 5% margin
        margin_x = 0.05 * (max(xs) - min(xs) or 1e-9)
        margin_y = 0.05 * (max(ys) - min(ys) or 1e-9)
        ax.set_xlim(min(xs) - margin_x, max(xs) + margin_x)
        ax.set_ylim(min(ys) - margin_y, max(ys) + margin_y)
        ax.grid(True, alpha=0.3)
        title = f"{model[0]} slice {axis_norm}"
        if args.x is not None and axis_norm != 'yz':
            title += f" | x={args.x}"
        if args.y is not None and axis_norm != 'xz':
            title += f" | y={args.y}"
        if args.z is not None and axis_norm != 'xy':
            title += f" | z={args.z}"
        ax.set_title(title)
        fig.tight_layout()
        if args.out:
            fig.savefig(args.out, dpi=200)
            print(f"Saved: {args.out}")
        else:
            plt.show()
        return

    if not args.plot_curve and args.unknown is None:
        print("Either provide --unknown for solving mode or use --plot-curve for plotting mode.", file=sys.stderr)
        sys.exit(2)

    unknown = args.unknown
    knowns = {"x": args.x, "y": args.y, "z": args.z}
    if unknown is not None and knowns[unknown] is not None:
        print(f"--unknown {unknown} but a value for {unknown} was provided. Remove it.", file=sys.stderr)
        sys.exit(2)
    provided = {k for k,v in knowns.items() if v is not None}
    if len(provided) != 2:
        print("Provide exactly two of --x, --y, --z.", file=sys.stderr)
        sys.exit(2)

    try:
        if model[0] == "linear":
            a,b,c = map(float, model[1])
            sols = _solve_linear(unknown, a,b,c, x=args.x, y=args.y, z=args.z)
        else:
            coef = list(map(float, model[1]))
            sols = _solve_quadratic(unknown, coef, x=args.x, y=args.y, z=args.z)
    except ZeroDivisionError as e:
        print(f"Solve failure: {e}", file=sys.stderr)
        sys.exit(3)

    if not sols:
        print("No real solution.", file=sys.stderr)
        sys.exit(1)

    # Print solutions, one per line: "<var>=<value>"
    for s in sols:
        print(f"{unknown}={s:.12g}")

if __name__ == "__main__":
    main()