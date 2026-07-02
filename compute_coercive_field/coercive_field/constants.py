"""Unit conversion constants.

这里集中放单位换算常数，避免在计算公式中散落魔法数字。
"""

# 1 (e * Angstrom) * 1 (MV / cm) = 0.01 eV。
# 外场修正使用 F = U - E * Delta_mu_parallel，因此 E 用 MV/cm、
# Delta_mu_parallel 用 e*Angstrom 时，乘这个因子得到 eV。
EV_PER_EANG_MV_CM = 0.01

# 1 (e * Angstrom / Angstrom^3) = 1602.176634 uC/cm^2.
# Formula:
#   P [uC/cm^2] = UC_CM2_PER_EANG_A3 * mu [e*Angstrom] / Omega [Angstrom^3]
# Symbols:
#   P: polarization.
#   mu: dipole moment.
#   Omega: cell volume.
UC_CM2_PER_EANG_A3 = 1602.176634
