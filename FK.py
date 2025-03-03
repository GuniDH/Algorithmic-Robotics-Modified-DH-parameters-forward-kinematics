import sympy as sp                    
from math import pi   

class viper300:
    """ 
    Describe the Viper300 6DOF robotic arm and perform FK on it (with normal transformation matrices or with modified dh parameters)
    In our physical model, the arm is pointing against the y-axis, 
    z is upwards, and x is pointing at the table's right side
    """
         
    # MDH determines if we use modified dh parameters to calculate the transformation matrices
    def __init__ (self, config, lengths, MDH=False):
        self.config = config
        self.lengths = lengths
        self.MDH = MDH
        self.T = None
    
    def mdh_transform(self, alpha_i_minus_1, a_i_minus_1,  d_i, phi_i):
        return sp.Matrix([
            [   sp.cos(phi_i),                         -sp.sin(phi_i),                                   0,                     a_i_minus_1                   ],
            [   sp.sin(phi_i)*sp.cos(alpha_i_minus_1), sp.cos(phi_i)*sp.cos(alpha_i_minus_1),  -sp.sin(alpha_i_minus_1),    -d_i*sp.sin(alpha_i_minus_1)      ],
            [   sp.sin(phi_i)*sp.sin(alpha_i_minus_1), sp.cos(phi_i)*sp.sin(alpha_i_minus_1),  sp.cos(alpha_i_minus_1),    d_i*sp.cos(alpha_i_minus_1)      ],
            [           0,                                        0,                                       0,                              1                      ]
        ])
    
    def calculate_Tx(self):
        q0, q1, q2, q3, q4 = self.config
        l1, l2, l3, l4, l5, l6 = self.lengths
        T = None
        
        if not self.MDH:          
            T01 = sp.Matrix([[sp.cos(q0), -sp.sin(q0), 0, 0],
                             [sp.sin(q0),  sp.cos(q0), 0, 0],
                             [0,           0         , 1, l1],
                             [0,           0         , 0, 1 ]])
            
            TRANSLATION12 =sp.Matrix([[1, 0, 0, 0],
                                      [0, 1, 0, 0],
                                      [0, 0, 1, l2],
                                      [0, 0, 0, 1]])
            
            ROTATION12 = sp.Matrix([[sp.cos(q1), 0,   -sp.sin(q1),  0],
                                    [0,          1,       0,        0],
                                    [sp.sin(q1), 0,   sp.cos(q1),   0],
                                    [0,          0,       0,        1]])
            
            T12 = ROTATION12*TRANSLATION12
            
            T23 =  sp.Matrix([[sp.cos(q2), 0,   -sp.sin(q2), -l3],
                              [0,          1,       0,        0],
                              [sp.sin(q2), 0,   sp.cos(q2),   0],
                              [0,          0,       0,        1]])
            
            T34 = sp.Matrix([[sp.cos(q3),     sp.sin(q3), 0,  0],
                             [-sp.sin(q3),    sp.cos(q3), 0,  0],
                             [0,                 0,       1,  l4+l5],
                             [0,                 0,       0,  1]])
            
            TRANSLATION45 = sp.Matrix([[1, 0, 0, 0],
                                       [0, 1, 0, -l6],
                                       [0, 0, 1, 0],
                                       [0, 0, 0, 1]])
            
            ROTATION45=sp.Matrix([[1,      0,             0,             0],
                                  [0, sp.cos(q4),      sp.sin(q4),       0],
                                  [0, -sp.sin(q4),     sp.cos(q4),       0],
                                  [0,      0,             0,             1]]) 
                
            T45=ROTATION45*TRANSLATION45
            
            T = T01*T12*T23*T34*T45
            
        else:
            MDH_params = [
    #  alpha i-1(s) a i-1(s)               d_i(s) phi_i(s)
            (0,      0,                      l1,    q0),
            (pi/2,   0,                      0,     q1 + (pi/2 + sp.atan(l3/l2))),
            (0,      sp.sqrt(l2**2+l3**2),   0,     q2 - (pi/2 + sp.atan(l3/l2))),
            (pi/2,   0,                  -(l4+l5),  q3 + pi/2),
            (-pi/2,  0,                      0,     q4),
            (0,      l6,                     0,     0)
            ]
            
            T = sp.eye(4)
            for i in range(6):
                T_i = self.mdh_transform(*MDH_params[i])
                T = T * T_i
            
        self.T = T
    
    def get_xyz(self):
        # Calculating forward kinematics to get ee location
        self.calculate_Tx()
        x = sp.Matrix([0, 0, 0, 1])
        Tx = self.T * x
        return Tx


lengths =  [120  * 1e-3,
            300  * 1e-3,
            60   * 1e-3,
            200  * 1e-3,
            100  * 1e-3,
            200  * 1e-3]
robot_configuration=[pi/3,pi/4,pi/5,pi/6,pi/7]

arm_normal = viper300(config=robot_configuration, lengths=lengths)
arm_mdh    = viper300(config=robot_configuration, lengths=lengths, MDH=True)
ee_xyz_normal, ee_xyz_mdh = arm_normal.get_xyz(), arm_mdh.get_xyz()

print(f'EE xyz of normal FK: {ee_xyz_normal[0]} {ee_xyz_normal[1]} {ee_xyz_normal[2]}')
print(f'EE xyz of MDH FK:    {ee_xyz_mdh[0]} {ee_xyz_mdh[1]} {ee_xyz_mdh[2]}')
