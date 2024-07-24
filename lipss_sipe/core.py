"""
Core functions for the calculation of the efficacy factor of LIPSS formation.

Equation numbers refer to the following paper:
Bonse, J., Munz, M., & Sturm, H. (2005). Structure formation on the surface of indium phosphide irradiated by femtosecond laser pulses. Journal of Applied Physics, 97(1), 013538. https://doi.org/10.1063/1.1827919

The code is published along with the following paper:
Kaczmarek, D., Albert, T. J., Bonse, J., Erratum: "Structure formation on the surface of indium phosphide irradiated by femtosecond laser pulses". https://doi.org/10.1063/5.0222903
"""

import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib.gridspec as gridspec
import matplotlib as mpl

FIGSIZE = (10, 7)
FONTSIZE = 25

# plot style settings
mpl.rc("image", interpolation="none")
mpl.rc("xtick", labelsize=18)
mpl.rc("ytick", labelsize=18)
mpl.rc("axes", labelsize=20)
mpl.rc("axes", titlesize=18)
mpl.rc("lines", linewidth=5)
mpl.rc("legend", fontsize=20)

mpl.rcParams["axes.linewidth"] = 3
mpl.rcParams["xtick.major.width"] = 3
mpl.rcParams["ytick.major.width"] = 3
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["axes.titleweight"] = "bold"
mpl.rcParams["figure.titleweight"] = "bold"
mpl.rcParams["figure.labelweight"] = "bold"
mpl.rcParams["font.weight"] = "bold"


def eta(polarization, theta_laser, kappa_x, kappa_y, epsilon, f, s):
    """
    Efficacy Factor.
    Eq. (1) in the paper.

    :param polarization: str
        Polarization of the incident light. Either "ppol" or "spol".
    :param theta_laser: float
        Angle of incidence of the laser beam in deg.
    :param kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    :param kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    :param epsilon: float
        Complex dielectric function of the material.
    :param f: float
        Filling factor.
    :param s: float
        Shape factor.
    :return: float
        Efficacy factor.
    """
    return (
        2
        * np.pi
        * np.absolute(
            v(polarization, theta_laser, kappa_x, kappa_y, epsilon, f, s, sign="p")
            + np.conj(
                v(polarization, theta_laser, kappa_x, kappa_y, epsilon, f, s, sign="m")
            )
        )
    )


def v(polarization, theta_laser, kappa_x, kappa_y, epsilon, f, s, sign):
    """
    Complex function v.
    Eq. (2a) for s-polarization and eq. (2b) for p-polarization in the paper.

    param: polarization: str
        Polarization of the incident light. Either "ppol" or "spol".
    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    param: kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    param: f: float
        Filling factor.
    param: s: float
        Shape factor.
    param: sign: str
        Sign of kappa_pm. Either "p" or "m".
    return: np.complex
        Complex function v.
    """
    if polarization == "ppol":
        return (
            (
                h_ss(theta_laser, kappa_x, kappa_y, epsilon, sign)
                * inner_prod_kx(theta_laser, kappa_x, kappa_y, sign) ** 2
                + h_kk(theta_laser, kappa_x, kappa_y, epsilon, sign)
                * inner_prod_ky(theta_laser, kappa_x, kappa_y, sign) ** 2
            )
            * gamma_t(epsilon, f, s)
            * np.absolute(t_x(theta_laser, epsilon)) ** 2
            + h_kz(theta_laser, kappa_x, kappa_y, epsilon, sign)
            * inner_prod_ky(theta_laser, kappa_x, kappa_y, sign)
            * gamma_z(epsilon, f, s)
            * epsilon
            * np.conj(t_x(theta_laser, epsilon))
            * t_z(theta_laser, epsilon)
            + h_zk(theta_laser, kappa_x, kappa_y, epsilon, sign)
            * inner_prod_ky(theta_laser, kappa_x, kappa_y, sign)
            * gamma_t(epsilon, f, s)
            * t_x(theta_laser, epsilon)
            * np.conj(t_z(theta_laser, epsilon))
            + h_zz(theta_laser, kappa_x, kappa_y, epsilon, sign)
            * gamma_z(epsilon, f, s)
            * epsilon
            * np.absolute(t_z(theta_laser, epsilon)) ** 2
        )
    if polarization == "spol":
        return (
            (
                h_ss(theta_laser, kappa_x, kappa_y, epsilon, sign)
                * inner_prod_ky(theta_laser, kappa_x, kappa_y, sign) ** 2
                + h_kk(theta_laser, kappa_x, kappa_y, epsilon, sign)
                * inner_prod_kx(theta_laser, kappa_x, kappa_y, sign) ** 2
            )
            * gamma_t(epsilon, f, s)
            * np.absolute(t_s(theta_laser, epsilon)) ** 2
        )


def inner_prod_ky(theta_laser, kappa_x, kappa_y, sign):
    """
    Inner product of k_+- and y.
    Eq. (3) in the paper.

    :param theta_laser: float
        Angle of incidence of the laser beam in deg.
    :param kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    :param kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    :param sign: str
        Sign of kappa_pm. Either "p" or "m".
    :return: float
        Inner product of k_+- and y.
    """
    if sign == "p":
        return (np.sin(theta_laser) + kappa_y) / kappa_pm(
            theta_laser, kappa_x, kappa_y, sign
        )
    if sign == "m":
        return (np.sin(theta_laser) - kappa_y) / kappa_pm(
            theta_laser, kappa_x, kappa_y, sign
        )
    else:
        raise ValueError("Invalid sign.")


def inner_prod_kx(theta_laser, kappa_x, kappa_y, sign):
    """
    Inner product of k_+- and x.
    Eq. (4) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    param: kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    param: sign: str
        Sign of kappa_pm. Either "p" or "m".
    return: float
        Inner product of k_+- and x.
    """
    if sign == "p":
        return +kappa_x / kappa_pm(theta_laser, kappa_x, kappa_y, sign)
    if sign == "m":
        return -kappa_x / kappa_pm(theta_laser, kappa_x, kappa_y, sign)
    else:
        raise ValueError("Invalid sign.")


def kappa(k, lambda_laser):
    """
    Dimensionless LIPSS wavevector.

    param: k: float
        Wavevector of the incident light.
    param: lambda_laser: float
        Wavelength of the incident light, in nm.
    return: float
        Dimensionless LIPSS wavevector.
    """
    return k * lambda_laser / (2 * np.pi)


def kappa_pm(theta_laser, kappa_x, kappa_y, sign):
    """
    Definition of kappa_pm.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    param: kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    param: sign: str
        Sign of kappa_pm. Either "p" or "m".
    return: float
        Definition of kappa_pm.
    """
    if sign == "p":
        return np.emath.sqrt(kappa_x**2 + (np.sin(theta_laser) + kappa_y) ** 2)
    if sign == "m":
        return np.emath.sqrt(kappa_x**2 + (np.sin(theta_laser) - kappa_y) ** 2)
    else:
        raise ValueError("Invalid sign.")


def h_ss(theta_laser, kappa_x, kappa_y, epsilon, sign):
    """
    Function h_ss.
    Eq. (5) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    param: kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    param: sign: str
        Sign of kappa_pm. Either "p" or "m".
    return: np.complex
        Function h_ss.
    """
    return 2j / (
        np.emath.sqrt(1 - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
        + np.emath.sqrt(epsilon - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
    )


def h_kk(theta_laser, kappa_x, kappa_y, epsilon, sign):
    """
    Function h_kk.
    Eq. (6) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    param: kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    param: sign: str
        Sign of kappa_pm. Either "p" or "m".
    return: np.complex
        Function h_kk.
    """
    return (
        2j
        * np.emath.sqrt((epsilon - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2))
        * np.emath.sqrt((1 - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2))
    ) / (
        epsilon * np.emath.sqrt(1 - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
        + np.emath.sqrt(epsilon - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
    )


def h_kz(theta_laser, kappa_x, kappa_y, epsilon, sign):
    """
    Function h_kz.
    Eq. (7) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    param: kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    param: sign: str
        Sign of kappa_pm. Either "p" or "m".
    return: np.complex
        Function h_kz.
    """
    return (
        2j
        * kappa_pm(theta_laser, kappa_x, kappa_y, sign)
        * np.emath.sqrt(epsilon - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
    ) / (
        epsilon * np.emath.sqrt(1 - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
        + np.emath.sqrt(epsilon - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
    )


def h_zk(theta_laser, kappa_x, kappa_y, epsilon, sign):
    """
    Function h_zk.
    Eq. (8) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    param: kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    param: sign: str
        Sign of kappa_pm. Either "p" or "m".
    return: np.complex
        Function h_zk.
    """
    return (
        2j
        * kappa_pm(theta_laser, kappa_x, kappa_y, sign)
        * np.emath.sqrt(1 - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
    ) / (
        epsilon * np.emath.sqrt(1 - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
        + np.emath.sqrt(epsilon - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
    )


def h_zz(theta_laser, kappa_x, kappa_y, epsilon, sign):
    """
    Function h_zz.
    Eq. (9) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: kappa_x: float
        Normalized LIPSS wavevector component in x-direction.
    param: kappa_y: float
        Normalized LIPSS wavevector component in y-direction.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    param: sign: str
        Sign of kappa_pm. Either "p" or "m".
    return: np.complex
        Function h_zz.
    """
    return (2j * kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2) / (
        epsilon * np.emath.sqrt(1 - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
        + np.emath.sqrt(epsilon - kappa_pm(theta_laser, kappa_x, kappa_y, sign) ** 2)
    )


def t_s(theta_laser, epsilon):
    """
    Complex function t_s.
    Eq. (10) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    return: np.complex
        Complex function t_s.
    """
    return (2 * np.absolute(np.cos(theta_laser))) / (
        np.abs(np.cos(theta_laser)) + np.emath.sqrt(epsilon - np.sin(theta_laser) ** 2)
    )


def t_x(theta_laser, epsilon):
    """
    Complex function t_x.
    Eq. (11) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    return: np.complex
        Complex function t_x.
    """
    return (
        2
        * np.emath.sqrt(epsilon - np.sin(theta_laser) ** 2)
        * np.abs(np.cos(theta_laser))
    ) / (
        epsilon * np.absolute(np.cos(theta_laser))
        + np.emath.sqrt(epsilon - np.sin(theta_laser) ** 2)
    )


def t_z(theta_laser, epsilon):
    """
    Complex function t_z.
    Eq. (12) in the paper.

    param: theta_laser: float
        Angle of incidence of the laser beam in deg.
    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    return: np.complex
        Complex function t_z.
    """
    return (2 * np.sin(theta_laser) * np.absolute(np.cos(theta_laser))) / (
        epsilon * np.absolute(np.cos(theta_laser))
        + np.emath.sqrt(epsilon - np.sin(theta_laser) ** 2)
    )


def gamma_t(epsilon, f, s):
    """
    Factor gamma_t which includes the surface roughness together with gamma_z.

    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    param: f: float
        Filling factor.
    param: s: float
        Shape factor.
    return: np.complex
        Factor gamma_t.
    """
    return (epsilon - 1) / (
        4 * np.pi * (1 + 0.5 * (1 - f) * (epsilon - 1) * (F(s) - R(epsilon) * G(s)))
    )


def gamma_z(epsilon, f, s):
    """
    Factor gamma_z which includes the surface roughness together with gamma_t.

    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    param: f: float
        Filling factor.
    param: s: float
        Shape factor.
    return: np.complex
        Factor gamma_z.
    """
    return (epsilon - 1) / (
        4 * np.pi * (epsilon - (1 - f) * (epsilon - 1) * (F(s) + R(epsilon) * G(s)))
    )


def R(epsilon):
    """
    Reflection coefficient R(epsilon).

    param: epsilon: np.complex
        Complex dielectric function of the material at the irradiation wavelength.
    return: np.complex
        Reflection coefficient R(epsilon).
    """
    return (epsilon - 1) / (epsilon + 1)


def F(s):
    """
    Scalar function F(s).
    Eq. (13) in the paper.

    param: s: float
        Shape factor.
    return: float
        Scalar function F(s).
    """
    return np.emath.sqrt(s**2 + 1) - s


def G(s):
    """
    Scalar function G(s).
    Eq. (14) in the paper.

    param: s: float
        Shape factor.
    return: float
        Scalar function G(s).
    """
    return 0.5 * (np.emath.sqrt(s**2 + 4) + s) - np.emath.sqrt(s**2 + 1)


def compute_eta(args):
    """
    Compute the efficacy factor eta for a single pair of wavevector components kappa_x and kappa_y.

    :param: args: tuple
        Tuple of arguments for the compute_eta function, defined as follows:
        i: int
            Index of the x-component of the wavevector
        k: int
            Index of the y-component of the wavevector
        x: float
            x-component of the wavevector
        y: float
            y-component of the wavevector
        polarization: str
            Polarization of the laser, either "spol" or "ppol"
        theta_laser: float
            Angle of incidence of the laser in radians
        epsilon: complex
            Complex dielectric function of the material at the irradiation wavelength
        f: float
            Filling factor
        s: float
            Shape factor
    :return: tuple
        Tuple of the input arguments and the computed efficacy factor eta, defined as follows:
        i: int
            Index of the x-component of the wavevector
        k: int
            Index of the y-component of the wavevector
        eta: float
            Efficacy factor eta
    """
    i, k, x, y, polarization, theta_laser, epsilon, f, s = args
    return (
        i,
        k,
        eta(
            polarization=polarization,
            theta_laser=theta_laser,
            kappa_x=x,
            kappa_y=y,
            epsilon=epsilon,
            f=f,
            s=s,
        ),
    )


def compute_eta_array(
    KAPPA_X, KAPPA_Y, POLARIZATION, THETA_LASER, EPSILON, F_FACTOR, S_FACTOR
):
    """
    Compute the efficacy factor eta for a range of wavevector components kappa_x and kappa_y.

    :param: KAPPA_X: np.ndarray
        1D array of the x-component of the wavevector
    :param: KAPPA_Y: np.ndarray
        1D array of the y-component of the wavevector
    :param: POLARIZATION: str
        Polarization of the laser, either "spol" or "ppol"
    :param: THETA_LASER: float
        Angle of incidence of the laser in radians
    :param: EPSILON: complex
        Complex dielectric function of the material at the irradiation wavelength
    :param: F_FACTOR: float
        Filling factor
    :param: S_FACTOR: float
        Shape factor
    :return: np.ndarray
        2D array of the efficacy factor eta
    """
    ETA = np.zeros((len(KAPPA_X), len(KAPPA_Y)), dtype=np.float64)

    # Create a list of arguments for the compute_eta function
    args_mp = [
        (i, k, x, y, POLARIZATION, THETA_LASER, EPSILON, F_FACTOR, S_FACTOR)
        for i, x in enumerate(KAPPA_X)
        for k, y in enumerate(KAPPA_Y)
    ]

    # Create a pool of worker processes
    with mp.Pool(mp.cpu_count()) as pool:
        # Use the pool to compute the eta values in parallel
        results = list(tqdm(pool.imap(compute_eta, args_mp), total=len(args_mp)))

    # Assign the computed eta values to the ETA array
    for i, k, eta_value in results:
        ETA[i, k] = eta_value

    return ETA


def plot_eta(ETA, KAPPA_X, KAPPA_Y):
    """
    Plot the efficacy factor eta as a function of the wavevector components kappa_x and kappa_y.

    :param: ETA: np.ndarray
        2D array of the efficacy factor eta
    :param: KAPPA_X: np.ndarray
        1D array of the x-component of the wavevector
    :param: KAPPA_Y: np.ndarray
        1D array of the y-component of the wavevector
    :return: None
    """
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    f = plt.figure(figsize=(12, 15))

    vmin = np.nanpercentile(ETA, 1)
    vmax = np.nanpercentile(ETA, 98)

    ax0 = plt.subplot(gs[0])

    im1 = ax0.pcolormesh(
        KAPPA_X,
        KAPPA_Y,
        ETA,
        vmin=vmin,
        vmax=vmax,
        cmap="plasma",
        shading="gouraud",
        rasterized=True,
        antialiased=True,
    )

    if len(KAPPA_X) == len(KAPPA_Y):
        ax0.set_aspect("equal")

    ax0.set_xlabel("$\kappa_x$", fontsize=20)
    ax0.set_ylabel("$\kappa_y$", fontsize=20)

    cbar = f.colorbar(im1, ax=ax0)
    cbar.set_label("Efficacy factor $\eta$", fontsize=20)

    ax1 = plt.subplot(gs[1])

    ax1.plot(
        KAPPA_X[np.where(KAPPA_X == 0)[0][0] :],
        ETA[np.where(KAPPA_X == 0)[0][0] :, np.where(KAPPA_X == 0)[0][0]],
        label=f"$\psi=0°$",
    )
    if is_square(ETA):
        zero_x_index = np.where(KAPPA_X == 0)[0][0]
        diagonal = np.diag(ETA)
        scale_factor = np.sqrt(1 + np.tan((np.pi / 2) - np.deg2rad(45)))
        ax1.plot(
            KAPPA_X[:] * scale_factor,
            diagonal,
            label=f"$\psi=45°$",
        )
    ax1.plot(
        KAPPA_Y[np.where(KAPPA_Y == 0)[0][0] :],
        ETA[np.where(KAPPA_Y == 0)[0][0], np.where(KAPPA_Y == 0)[0][0] :],
        label=f"$\psi=90°$",
    )
    ax1.set_xlabel("$\kappa$", fontsize=20)
    ax1.set_ylabel("$\eta$", fontsize=20)
    ax1.set_xlim(0, max(KAPPA_X))
    ax1.legend()

    plt.tight_layout()
    plt.show()


def is_square(matrix):
    """
    Check if a matrix is square.

    :param: matrix: np.ndarray
        2D array
    :return: bool
        True if the matrix is square, False otherwise
    """
    return matrix.shape[0] == matrix.shape[1]
