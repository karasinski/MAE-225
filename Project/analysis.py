from sympy import *
from sympy.matrices import *

class Transform:
    def __init__(self, phi, alpha, a, d):
        self.phi = phi
        self.alpha = alpha
        self.a = a
        self.d = d
        self.DHMatrix()

    def DHMatrix(self):
        phi = rad(self.phi)
        alpha = rad(self.alpha)
        a = self.a
        d = self.d

        dh = Matrix([[cos(phi), -sin(phi)*cos(alpha)  ,sin(phi)*sin(alpha), a*cos(phi)],
                     [sin(phi),  cos(phi)*cos(alpha), -cos(phi)*sin(alpha), a*sin(phi)],
                     [       0,           sin(alpha),           cos(alpha),          d],
                     [       0,                    0,                    0,          1]])
        self.A = dh

        self.R = dh[:3, :3]
        self.vecA = dh[:3, 3]

# Let's do an analysis with the SSRMS with the shoulder roll joint locked at an angle beta

beta = symbols('beta')
Theta1, Theta2, Theta3, Theta4, Theta5, Theta6, Theta7 = symbols("Theta1 Theta2 Theta3 Theta4 Theta5 Theta6 Theta7")
d10, d20, d30, d60, d70, a30, a40 = symbols('d10 d20 d30 d60 d70 a30 a40')

# # # T1 = DHMatrix( 90 + Theta1, 90,   0, d10);
# # T2 = DHMatrix( 90 + Theta2 + beta, 90,   0, d20);
# # T3 = DHMatrix(  0 + Theta3,         0, a30, d30);
# # T4 = DHMatrix(  0 + Theta4,         0, a40,   0);
# # T5 = DHMatrix(180 + Theta5,        90,   0,   0);
# # T6 = DHMatrix(-90 + Theta6,        90,   0, d60);
# # T7 = DHMatrix(180 + Theta7,        90,   0, d70);

# # Define the DH matrices for each link
# # T1 = DHMatrix(Theta1, 90,   0, d10);
# A23 = DHMatrix(Theta2 + beta, 90,   0, d20);
# A34 = DHMatrix(Theta3,         0, a30, d30);
# A45 = DHMatrix(Theta4,         0, a40,   0);
# A56 = DHMatrix(Theta5,        90,   0,   0);
# A67 = DHMatrix(Theta6,        90,   0, d60);
# A78 = DHMatrix(Theta7,        90,   0, d70);

# # Grab the rotation matrix from each A matrix
# R23 = A23[:3, :3]
# R34 = A34[:3, :3]
# R45 = A45[:3, :3]
# R56 = A56[:3, :3]
# R67 = A67[:3, :3]
# R78 = A78[:3, :3]

# # And grab the displacement vectors
# a2 = A23[:3, 3]
# a3 = A34[:3, 3]
# a4 = A45[:3, 3]
# a5 = A56[:3, 3]
# a6 = A67[:3, 3]
# a7 = A78[:3, 3]

# # The z_i's need to be WRT the base coordinate system
# z = Matrix([[0], [0], [1]])

# z2 = z
# z3 = R23 * z
# z4 = R23 * R34 * z
# z5 = R23 * R34 * R45 * z
# z6 = R23 * R34 * R45 * R56 * z
# z7 = R23 * R34 * R45 * R56 * R67 * z

# r7 = R67 * 

# # Let's assume that the shoulder roll joint is locked at beta = 0
# T = (A2 * A3 * A4 * A5 * A6 * A7).subs({beta: 0})


# configs = [[Theta2 + beta, 90,   0, d20],
#          [Theta3,         0, a30, d30],
#          [Theta4,         0, a40,   0],
#          [Theta5,        90,   0,   0],
#          [Theta6,        90,   0, d60],
#          [Theta7,        90,   0, d70]]


configs = [[Theta1, -90, 0,   0],
           [Theta2,  90, 0, d20],
           [   -90,   0, 0, d30],
           [Theta4, -90, 0,   0],
           [Theta5,  90, 0,   0],
           [Theta6,   0, 0,   0]]

T = []
for config in configs:
    T.append(Transform(*config))

# T23 = Transform(Theta2 + beta, 90,   0, d20);
# T34 = Transform(Theta3,         0, a30, d30);
# T45 = Transform(Theta4,         0, a40,   0);
# T56 = Transform(Theta5,        90,   0,   0);
# T67 = Transform(Theta6,        90,   0, d60);
# T78 = Transform(Theta7,        90,   0, d70);

# The z_i's need to be WRT the base coordinate system
z = Matrix([[0], [0], [1]])

for i, _ in enumerate(T):
    rot = Identity(3)
    for j in range(i):
        rot *= T[j].R 
    T[i].z = rot * z

# T23.z = z
# T34.z = T23.R * z
# T45.z = T23.R * T34.R * z
# T56.z = T23.R * T34.R * T45.R * z
# T67.z = T23.R * T34.R * T45.R * T56.R * z
# T78.z = T23.R * T34.R * T45.R * T56.R * T67.R * z

for i, _ in enumerate(T):
    i = len(T) - (i+1)
    rot = Identity(3)
    for j in range(i):
        rot *= T[j].R
    rot *= T[i].vecA
    if (i+1) < len(T):
        rot += T[i+1].r
    T[i].r = rot

# T78.r = T23.R * T34.R * T45.R * T56.R * T67.R * T78.vecA
# T67.r = T23.R * T34.R * T45.R * T56.R * T67.vecA + T78.r
# T56.r = T23.R * T34.R * T45.R * T56.vecA + T67.r
# T45.r = T23.R * T34.R * T45.vecA + T56.r
# T34.r = T23.R * T34.vecA + T45.r
# T23.r = T23.vecA + T34.r

for i, _ in enumerate(T):
    T[i].z_skew = Matrix([[         0, -T[i].z[2],  T[i].z[1]],
                          [ T[i].z[2],          0, -T[i].z[0]],
                          [-T[i].z[1],  T[i].z[0],         0]])
    T[i].crossed = T[i].z_skew * T[i].r

top, bottom = Matrix(), Matrix()
for i, _ in enumerate(T):
    top = top.row_join(T[i].crossed)
    bottom = bottom.row_join(T[i].z)
J = top.col_join(bottom)
