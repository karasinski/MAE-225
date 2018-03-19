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
beta, theta1, theta2, theta3, theta4, theta5, theta6, theta7 = symbols("beta theta1 theta2 theta3 theta4 theta5 theta6 theta7")
d1, d2, d3, d6, d7, a3, a4 = symbols('d1 d2 d3 d6 d7 a3 a4')

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

def KJ():
  # Kinematic Jacobian
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

  return J

def GJ():
  # Geometric Jacobian
  Q_rev = Matrix([[0, -1, 0, 0],
                  [1,  0, 0, 0],
                  [0,  0, 0, 0],
                  [0,  0, 0, 0]])

  for i, _ in enumerate(T):
    print('i ', i)
    if i > 0:
      T[i].Ti = T[i-1].Ti * T[i-1].A
    if i == 0:
      T[i].Ti = Matrix(Identity(4))
    T[i].Ti.simplify()

    T[i].D = simplify(T[i].Ti * Q_rev) * simplify(T[i].Ti.inv())

  Jg = Matrix([[T[0].D[0, 3], T[1].D[0, 3], T[2].D[0, 3], T[3].D[0, 3], T[4].D[0, 3], T[5].D[0, 3], T[6].D[0, 3]],
               [T[0].D[1, 3], T[1].D[1, 3], T[2].D[1, 3], T[3].D[1, 3], T[4].D[1, 3], T[5].D[1, 3], T[6].D[1, 3]],
               [T[0].D[2, 3], T[1].D[2, 3], T[2].D[2, 3], T[3].D[2, 3], T[4].D[2, 3], T[5].D[2, 3], T[6].D[2, 3]],
               [T[0].D[2, 1], T[1].D[2, 1], T[2].D[2, 1], T[3].D[2, 1], T[4].D[2, 1], T[5].D[2, 1], T[6].D[2, 1]],
               [T[0].D[0, 2], T[1].D[0, 2], T[2].D[0, 2], T[3].D[0, 2], T[4].D[0, 2], T[5].D[0, 2], T[6].D[0, 2]],
               [T[0].D[1, 0], T[1].D[1, 0], T[2].D[1, 0], T[3].D[1, 0], T[4].D[1, 0], T[5].D[1, 0], T[6].D[1, 0]]])
  return Jg

# Inverse Kinematic Solution
# T07 = T[0].A * T[1].A * T[2].A * T[3].A * T[4].A * T[5].A * T[6].A

# We know that
# u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, p_x, p_y, p_z = symbols('u_x u_y u_z v_x v_y v_z w_x w_y w_z p_x p_y p_z')
# T07 = Matrix([[u_x, v_x, w_x, p_x],
#               [u_y, v_y, w_y, p_y],
#               [u_z, v_z, w_z, p_z],
#               [  0,   0,   0,   1]])

# And setting T01^-1 * T07 = T12 T23 T34 T45 T56 T67
# print(latex(T[0].A.inv() * T07))
# print(latex(T[1].A * T[2].A * T[3].A * T[4].A * T[5].A * T[6].A))

def IK(SHOULDER=1, WRIST=1, ELBOW=1,
       a3=2.3, a4=2.3, d1=0.65, d2=0.3,
       d3=0.9, d6=0.3, d7=0.65,
       beta=rad(60),
       u_x=0.8021, u_y=-0.5859, u_z=-0.1154,
       v_x=0.1217, v_y= 0.3495, v_z=-0.9290,
       w_x=0.5846, w_y= 0.7311, w_z=0.3517,
       p_x=2.4790, p_y=-2.4734, p_z=-0.4927):
  '''
  Return the IK for a given pose, initial locked angle, and arm characteristics.
  '''

  theta1 = beta

  h1 = -v_z * d7 - d1 + p_z
  q1 = (v_x * d7 - p_x)*cos(beta) + (v_y * d7 - p_y)*sin(beta)

  theta2 = SHOULDER * acos(d3/(sqrt(h1**2 + q1**2))) + atan2(q1, h1) - pi

  theta6 = WRIST * acos(v_z*cos(theta2)-(v_x*cos(theta1) + v_y*sin(theta1))*sin(theta2))

  theta7 = -atan2(((u_x*cos(theta1) + u_y*sin(theta1))*sin(theta2) - u_z*cos(theta2))/sin(theta6),
                  ((w_x*cos(theta1) + w_y*sin(theta1))*sin(theta2) - w_z*cos(theta2))/sin(theta6)) + pi/2

  X = d6 * ((w_z * sin(theta2) + cos(theta2) * (w_x * cos(theta1) + w_y * sin(theta1))) * cos(theta7) - (u_z * sin(theta2) + cos(theta2) * (u_x * cos(theta1) + u_y * sin(theta1))) * sin(theta7)) - d7 * (v_z * sin(theta2) + cos(theta2) * (v_x * cos(theta1) + v_y * sin(theta1))) + (-d1 + p_z) * sin(theta2) + cos(theta2) * (p_x * cos(theta1) + p_y * sin(theta1))
  Y = -d2 + d6 * ((w_x * sin(theta1) - w_y * cos(theta1)) * cos(theta7) - (u_x * sin(theta1) - u_y * cos(theta1)) * sin(theta7)) - d7 * (v_x * sin(theta1) - v_y * cos(theta1)) + p_x * sin(theta1) - p_y*cos(theta1)

  theta4 = ELBOW * acos((X**2 + Y**2 - a3**2 - a4**2)/(2*a3*a4))

  theta3 = atan2(Y*(a3 + a4*cos(theta4)) - X*a4*sin(theta4),
                 X*(a3 + a4*cos(theta4)) + Y*a4*sin(theta4))

  theta5 = atan2((v_x*sin(theta1)-v_y*cos(theta1))/(sin(theta6)),
                 (w_x*sin(theta1)-w_y*cos(theta1))*cos(theta7)-(u_x*sin(theta1)-u_y*cos(theta1))*sin(theta7)) - (theta3 + theta4)

  thetas = []
  for theta in [theta1, theta2, theta3, theta4, theta5, theta6, theta7]:
    thetas.append(deg(theta).evalf())

  thetas = Matrix(thetas)
  return thetas

def IK_sols():
  '''
  Enumerate through the eight possible solutions for a given pose.
  '''

  c = [[ 1,  1,  1],
       [ 1,  1, -1],
       [ 1, -1,  1],
       [-1,  1,  1],
       [ 1, -1, -1],
       [-1,  1, -1],
       [-1, -1,  1],
       [-1, -1, -1]]

  res = []
  for S, E, W in c:
      print(S, E, W)
      t = IK(SHOULDER=S, ELBOW=E, WRIST=W)
      res.append(t)

  stacked = np.zeros((7, 1))
  for row in res:
      stacked = np.hstack((stacked, np.array(row)))

  return stacked[:, 1:].T
