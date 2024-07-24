import argparse
import numpy as np
from lipss_sipe.core import compute_eta_array, plot_eta


def main():
    """
    Command line interface for the LIPSS-SIPE package.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--f", type=float, help="Filling factor")
    parser.add_argument("--s", type=float, help="Shape factor")
    parser.add_argument(
        "--n",
        type=float,
        nargs=2,
        metavar=("n", "k"),
        help="Refractive index n and absorption coefficient k at the irradiation wavelength",
    )
    # parser.add_argument("--wavelength", type=float, help="Laser wavelength in nm")
    parser.add_argument(
        "--angle", type=float, help="Angle of incidence of the laser in degrees"
    )
    parser.add_argument(
        "--polarization",
        type=str,
        choices=("s", "p"),
        help="Polarization of the laser",
    )
    parser.add_argument(
        "--x_range",
        type=int,
        nargs=3,
        metavar=("start", "stop", "num"),
        help="Range and number of normalized values for the x-component of the wavevector",
    )
    parser.add_argument(
        "--y_range",
        type=int,
        nargs=3,
        metavar=("start", "stop", "num"),
        help="Range and number of normalized values for the y-component of the wavevector",
    )

    args = parser.parse_args()

    # if a command line argument is not provided, ask the user for input
    if args.f is None:
        args.f = float(input("Enter the filling factor (for example 0.1): "))
    if args.s is None:
        args.s = float(input("Enter the shape factor (for example 10): "))
    if args.n is None:
        args.n = input(
            "Enter the refractive index n and the absorption coefficient k (for example 4 0.1): "
        )
        args.n = [float(x) for x in args.n.split()]
    # if args.wavelength is None:
    #     args.wavelength = float(input("Enter the laser wavelength in nm (for example 1060): "))
    if args.angle is None:
        args.angle = float(
            input(
                "Enter the angle of incidence of the laser in degrees (for example 60): "
            )
        )
    if args.polarization is None:
        args.polarization = input(
            "Enter the polarization of the laser, either 's' or 'p': "
        )
    if args.x_range is None:
        args.x_range = input(
            "Enter the range and number of values for the x-component of the wavevector (for example -5 5 151): "
        )
        args.x_range = [int(x) for x in args.x_range.split()]
    if args.y_range is None:
        args.y_range = input(
            "Enter the range and number of values for the y-component of the wavevector (for example -5 5 151): "
        )
        args.y_range = [int(x) for x in args.y_range.split()]

    # polarization of the laser, either "spol" or "ppol"
    if args.polarization == "s":
        POLARIZATION = "spol"
    elif args.polarization == "p":
        POLARIZATION = "ppol"
    else:
        raise ValueError("Invalid polarization. Please enter 's' or 'p'.")

    # surface roughness properties
    F_FACTOR = args.f  # filling factor
    S_FACTOR = args.s  # shape factor

    EPSILON = (
        args.n[0] + 1j * args.n[1]
    ) ** 2  # complex dielectric function of the material at the irradiation wavelength
    # LAMBDA_LASER = args.wavelength  # laser wavelength
    THETA_LASER = np.deg2rad(args.angle)  # angle of incidence of the laser

    KAPPA_X = np.linspace(
        args.x_range[0], args.x_range[1], args.x_range[2]
    )  # x-component of the wavevector
    KAPPA_Y = np.linspace(
        args.y_range[0], args.y_range[1], args.x_range[2]
    )  # y-component of the wavevector

    ETA = compute_eta_array(
        KAPPA_X, KAPPA_Y, POLARIZATION, THETA_LASER, EPSILON, F_FACTOR, S_FACTOR
    )

    plot_eta(ETA, KAPPA_X, KAPPA_Y)
