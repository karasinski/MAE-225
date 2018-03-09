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
        # phi = rad(self.phi)
        phi = self.phi
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
beta, theta1, theta2, theta3, theta4, theta5, theta6, theta7 = symbols("beta theta1 theta2 theta3 theta4 theta5 theta6 theta7", real=True)
d1, d2, d3, d6, d7, a3, a4 = symbols('d1 d2 d3 d6 d7 a3 a4', positive=True)

# SSRMS
configs = [[theta1, 90,  0, d1],
           [theta2, 90,  0, d2],
           [theta3,  0, a3, d3],
           [theta4,  0, a4,  0],
           [theta5, 90,  0,  0],
           [theta6, 90,  0, d6],
           [theta7, 90,  0, d7]]

# Stanford Arm Example from Class
# configs = [[theta1, -90, 0,  0],
#            [theta2,  90, 0, d2],
#            [   -90,   0, 0, d3],
#            [theta4, -90, 0,  0],
#            [theta5,  90, 0,  0],
#            [theta6,   0, 0,  0]]

# Build our robot
T = []
for config in configs:
    T.append(Transform(*config))

# The z_i's need to be WRT the base coordinate system
z = Matrix([[0], [0], [1]])

# Calculate all our z vectors
for i, _ in enumerate(T):
    rot = Identity(3)
    for j in range(i):
        rot *= T[j].R 
    T[i].z = rot * z

# Loop through backwards to calculate the r vectors
for i, _ in enumerate(T):
    i = len(T) - (i+1)
    rot = Identity(3)
    for j in range(i):
        rot *= T[j].R
    rot *= T[i].vecA
    if (i+1) < len(T):
        rot += T[i+1].r
    T[i].r = rot

# Calculate the z skew matrices and compute cross products
for i, _ in enumerate(T):
    T[i].z_skew = Matrix([[         0, -T[i].z[2],  T[i].z[1]],
                          [ T[i].z[2],          0, -T[i].z[0]],
                          [-T[i].z[1],  T[i].z[0],         0]])
    T[i].crossed = T[i].z_skew * T[i].r

# Combine elements to form Jacobian
top, bottom = Matrix(), Matrix()
for i, _ in enumerate(T):
    top = top.row_join(T[i].crossed)
    bottom = bottom.row_join(T[i].z)
J = top.col_join(bottom)

# Assume that beta is set to zero
# J0 = J.subs({beta: 0})

# T07 = T[0].A * T[1].A * T[2].A * T[3].A * T[4].A * T[5].A * T[6].A

# We know that 
n_x, n_y, n_z, o_x, o_y, o_z, a_x, a_y, a_z, p_x, p_y, p_z = symbols('n_x n_y n_z o_x o_y o_z a_x a_y a_z p_x p_y p_z')
T07 = Matrix([[n_x, o_x, a_x, p_x],
              [n_y, o_y, a_y, p_y],
              [n_z, o_z, a_z, p_z],
              [  0,   0,   0,   1]])

# And setting T01^-1 * T07 = T12 T23 T34 T45 T56 T67
# print(latex(T[0].A.inv() * T07))

Q_rev = Matrix([[0, -1, 0, 0],
                [1,  0, 0, 0],
                [0,  0, 0, 0],
                [0,  0, 0, 0]])

for i, _ in enumerate(T):
  print('i ', i)
  if i > 0:
    T[i].Ti = T[i-1].Ti * T[i].A
  if i == 0:
    T[i].Ti = T[i].A
  T[i].Ti.simplify()

  T[i].D = simplify(T[i].Ti * Q_rev) * simplify(T[i].Ti.inv())
  # print(latex(Di))

for i in range(7):
  print('i ', i)
  print(latex(simplify(T[i].D[0, 3])))
  print(latex(simplify(T[i].D[1, 3])))
  print(latex(simplify(T[i].D[2, 3])))
  print(latex(simplify(T[i].D[2, 1])))
  print(latex(simplify(T[i].D[0, 2])))
  print(latex(simplify(T[i].D[1, 0])))
  print('\n')

# ('i ', 0)                 
# ('i ', 1)                 
# ('i ', 2)                 
# ('i ', 3)                 
# ('i ', 4)                 
# ('i ', 5)                 
# ('i ', 6)                 

# ('i ', 0)                                            
# d_{1} \cos{\left (\theta_{1} \right )}               
# d_{1} \sin{\left (\theta_{1} \right )}               
# 0                                                    
# \sin{\left (\theta_{1} \right )}                     
# - \cos{\left (\theta_{1} \right )}                   
# 0                                                    


# ('i ', 1)                                            
# - d_{1} \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )} + d_{2} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}                                                                  
# d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )} + d_{2} \sin{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}                                                                    
# d_{2} \sin{\left (\theta_{2} \right )}               
# \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )}                                         
# \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )}                                         
# - \cos{\left (\theta_{2} \right )}                   


# ('i ', 2)                                            
# - a_{3} \sin{\left (\theta_{1} \right )} \cos{\left (\theta_{3} \right )} + a_{3} \sin{\left (\theta_{3} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} - d_{1} \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )} + d_{2} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}
# a_{3} \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} \right )} \cos{\left (\theta_{2} \right )} + a_{3} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{3} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )} + d_{2} \sin{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}
# \left(a_{3} \sin{\left (\theta_{3} \right )} + d_{2}\right) \sin{\left (\theta_{2} \right )}              
# \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )}                                         
# \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )}                                         
# - \cos{\left (\theta_{2} \right )}                   


# ('i ', 3)                                            
# - a_{3} \sin{\left (\theta_{1} \right )} \cos{\left (\theta_{3} \right )} + a_{3} \sin{\left (\theta_{3} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} - a_{4} \sin{\left (\theta_{1} \right )} \cos{\left (\theta_{3} + \theta_{4} \right )} + a_{4} \sin{\left (\theta_{3} + \theta_{4} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} - d_{1} \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )} + d_{2} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}
# a_{3} \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} \right )} \cos{\left (\theta_{2} \right )} + a_{3} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{3} \right )} + a_{4} \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} \right )} \cos{\left (\theta_{2} \right )} + a_{4} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{3} + \theta_{4} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )} + d_{2} \sin{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}
# \left(a_{3} \sin{\left (\theta_{3} \right )} + a_{4} \sin{\left (\theta_{3} + \theta_{4} \right )} + d_{2}\right) \sin{\left (\theta_{2} \right )}                                                                   
# \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )}                                         
# \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )}                                         
# - \cos{\left (\theta_{2} \right )}                   


# ('i ', 4)                                            
# - \left(d_{1} \cos{\left (\theta_{2} \right )} - d_{3}\right) \left(\sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) - \left(a_{3} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \cos{\left (\theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )}
# - \left(d_{1} \cos{\left (\theta_{2} \right )} - d_{3}\right) \left(\sin{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )}\right) - \left(a_{3} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \cos{\left (\theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )}
# a_{3} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{5} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{2} \right )} + d_{3} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}
# - \sin{\left (\theta_{1} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}
# \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{2} \right )} + \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}
# \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}               


# ('i ', 5)                                            
# \left(\left(\sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) \cos{\left (\theta_{6} \right )} + \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{1} \right )}\right) \left(a_{3} \sin{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{2} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{6}\right) + \left(\sin{\left (\theta_{1} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}\right) \left(a_{3} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \cos{\left (\theta_{5} \right )} \cos{\left (\theta_{6} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{1} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{6} \right )} + d_{3} \sin{\left (\theta_{6} \right )}\right)
# \left(\left(\sin{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )}\right) \cos{\left (\theta_{6} \right )} + \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )}\right) \left(a_{3} \sin{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{2} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{6}\right) - \left(\sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{2} \right )} + \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) \left(a_{3} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \cos{\left (\theta_{5} \right )} \cos{\left (\theta_{6} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{1} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{6} \right )} + d_{3} \sin{\left (\theta_{6} \right )}\right)
# \left(\sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )}\right) \left(a_{3} \sin{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{2} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{6}\right) - \left(a_{3} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \cos{\left (\theta_{5} \right )} \cos{\left (\theta_{6} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{1} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{6} \right )} + d_{3} \sin{\left (\theta_{6} \right )}\right) \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}
# \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{6} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{6} \right )} + \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}
# - \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} + \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{6} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )}
# \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )}                     


# ('i ', 6)                                            
# \left(\left(\left(\sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) \cos{\left (\theta_{6} \right )} + \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{1} \right )}\right) \cos{\left (\theta_{7} \right )} - \left(\sin{\left (\theta_{1} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )}\right) \sin{\left (\theta_{7} \right )}\right) \left(a_{3} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{1} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} + d_{2} \sin{\left (\theta_{6} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{3} \cos{\left (\theta_{6} \right )} + d_{7}\right) - \left(\left(\sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) \sin{\left (\theta_{6} \right )} - \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{6} \right )}\right) \left(a_{3} \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{4} + \theta_{5} \right )} + a_{3} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{5} \right )} \sin{\left (\theta_{7} \right )} + a_{4} \cos{\left (\theta_{5} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{1} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{7} \right )} - d_{2} \sin{\left (\theta_{7} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{3} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{6} \sin{\left (\theta_{7} \right )}\right)
# \left(\left(\left(\sin{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )}\right) \cos{\left (\theta_{6} \right )} + \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )}\right) \cos{\left (\theta_{7} \right )} + \left(\sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{2} \right )} + \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) \sin{\left (\theta_{7} \right )}\right) \left(a_{3} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{1} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} + d_{2} \sin{\left (\theta_{6} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{3} \cos{\left (\theta_{6} \right )} + d_{7}\right) - \left(\left(\sin{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )}\right) \sin{\left (\theta_{6} \right )} - \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )}\right) \left(a_{3} \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{4} + \theta_{5} \right )} + a_{3} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{5} \right )} \sin{\left (\theta_{7} \right )} + a_{4} \cos{\left (\theta_{5} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{1} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{7} \right )} - d_{2} \sin{\left (\theta_{7} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{3} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{6} \sin{\left (\theta_{7} \right )}\right)
# \left(\left(\sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )}\right) \cos{\left (\theta_{7} \right )} + \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}\right) \left(a_{3} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{1} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} + d_{2} \sin{\left (\theta_{6} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{3} \cos{\left (\theta_{6} \right )} + d_{7}\right) - \left(\sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )}\right) \left(a_{3} \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{4} + \theta_{5} \right )} + a_{3} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} \cos{\left (\theta_{4} + \theta_{5} \right )} + a_{4} \sin{\left (\theta_{5} \right )} \sin{\left (\theta_{7} \right )} + a_{4} \cos{\left (\theta_{5} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{1} \sin{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - d_{1} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{7} \right )} - d_{2} \sin{\left (\theta_{7} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + d_{2} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{3} \sin{\left (\theta_{6} \right )} \cos{\left (\theta_{7} \right )} + d_{6} \sin{\left (\theta_{7} \right )}\right)
# \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{6} \right )} + \sin{\left (\theta_{1} \right )} \cos{\left (\theta_{7} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} + \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \sin{\left (\theta_{7} \right )} \cos{\left (\theta_{1} \right )} + \sin{\left (\theta_{7} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{7} \right )}
# \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{6} \right )} \sin{\left (\theta_{7} \right )} + \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{7} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{1} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{2} \right )} \cos{\left (\theta_{7} \right )} - \sin{\left (\theta_{7} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{6} \right )} - \cos{\left (\theta_{1} \right )} \cos{\left (\theta_{7} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )}
# \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{7} \right )} \cos{\left (\theta_{6} \right )} \cos{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} - \sin{\left (\theta_{2} \right )} \sin{\left (\theta_{3} + \theta_{4} + \theta_{5} \right )} \cos{\left (\theta_{7} \right )} - \sin{\left (\theta_{6} \right )} \sin{\left (\theta_{7} \right )} \cos{\left (\theta_{2} \right )}
