import IAPWS_IF97 as iapws

T = 647.1
P = 19.585

density = iapws.density(P, T)

print(density)
