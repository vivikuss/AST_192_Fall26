"""Step 8: Galactic UVW relative to the LSR and height z, with astropy."""
import astropy.units as u
import numpy as np
from astropy.coordinates import (ICRS, CartesianDifferential, Galactocentric,
                                 GalacticLSR, SkyCoord)


def skycoord(df, dist_col="dist_pc", rv_col="rv"):
    return SkyCoord(ra=np.asarray(df["ra"]) * u.deg, dec=np.asarray(df["dec"]) * u.deg,
                    distance=np.asarray(df[dist_col]) * u.pc,
                    pm_ra_cosdec=np.asarray(df["pmra"]) * u.mas / u.yr,
                    pm_dec=np.asarray(df["pmdec"]) * u.mas / u.yr,
                    radial_velocity=np.asarray(df[rv_col]) * u.km / u.s, frame=ICRS())


def uvw(c, solar_motion_kms=(11.1, 12.24, 7.25)):
    """U (to Galactic centre), V (rotation), W (NGP) in km/s relative to the LSR.

    solar_motion_kms = Sun's velocity w.r.t. the LSR (astropy GalacticLSR v_bary).
    """
    v_bary = CartesianDifferential(np.asarray(solar_motion_kms, dtype=float) * u.km / u.s)
    g = c.transform_to(GalacticLSR(v_bary=v_bary))
    v = g.velocity
    return v.d_x.to_value(u.km / u.s), v.d_y.to_value(u.km / u.s), v.d_z.to_value(u.km / u.s)


def uvw_heliocentric(c):
    """U, V, W relative to the Sun (no LSR correction), for literature checks."""
    return uvw(c, (0.0, 0.0, 0.0))


def height_z(c, kcfg):
    """Height above the Galactic mid-plane [pc]."""
    if kcfg["z_definition"] == "heliocentric":
        g = c.galactic
        return (g.distance * np.sin(g.b)).to_value(u.pc)
    frame = Galactocentric(galcen_distance=kcfg["galcen_distance_kpc"] * u.kpc,
                           z_sun=kcfg["z_sun_pc"] * u.pc)
    return c.transform_to(frame).z.to_value(u.pc)


def add_kinematics(df, cfg, dist_col="dist_pc", rv_col="rv"):
    c = skycoord(df, dist_col, rv_col)
    k = cfg["kinematics"]
    df["U"], df["V"], df["W"] = uvw(c, k["solar_motion_kms"])
    df["z"] = height_z(c, k)
    return df


def dispersion(U, V, W, statistic="rms"):
    """Velocity dispersion about the sample's own mean.

    rms      : sqrt( (1/N) sum[(U-Ubar)^2 + (V-Vbar)^2 + (W-Wbar)^2] )   (spec; D14)
    mean_abs : (1/N) sum[(U-Ubar)^2 + (V-Vbar)^2 + (W-Wbar)^2]^(1/2)   (HS19 eq. 1 as typeset)

    Accepts 1-D arrays, or 2-D (iterations x N) arrays -> one value per row.
    """
    U, V, W = (np.asarray(a, dtype=float) for a in (U, V, W))
    ax = -1
    ss = ((U - U.mean(ax, keepdims=True)) ** 2 + (V - V.mean(ax, keepdims=True)) ** 2
          + (W - W.mean(ax, keepdims=True)) ** 2)
    if statistic == "rms":
        return np.sqrt(ss.mean(ax))
    if statistic == "mean_abs":
        return np.sqrt(ss).mean(ax)
    raise ValueError(statistic)
