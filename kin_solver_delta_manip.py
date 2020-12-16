import math


class DeltaKinematics:
    def __init__(self, e, f, re, rf):
        """
        e - сторона подвижной платформы
        f - сторона неподвижного основания
        re - длина большего плеча
        rf - длина короткого плечв
        """
        # длины опор и звеньев
        self.e = e
        self.f = f
        self.re = re
        self.rf = rf

        # тригонометрические константы
        self.sqrt3 = math.sqrt(3.0)
        self.sin120 = self.sqrt3 / 2.0
        self.cos120 = -0.5
        self.tan60 = self.sqrt3
        self.sin30 = 0.5
        self.tan30 = 1 / self.sqrt3

        # координаты центра подвижного основания
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0

        # углы поворота звеньев
        self.theta_1 = 0.0
        self.theta_2 = 0.0
        self.theta_3 = 0.0

    def forward_kinematics(self, thetas):
        """
        прямая кинематика: (theta1, theta2, theta3) -> (x0, y0, z0)
        возвращаемый статус: 0=OK, -1=несуществующая позиция
        """
        t = (self.f - self.e) * self.tan30 / 2
        for i, theta in enumerate(thetas):
            thetas[i] = theta * math.pi / 180.0
    
        y_1 = -(t + self.rf * math.cos(thetas[0]))
        z_1 = -self.rf * math.sin(thetas[0])
    
        y_2 = (t + self.rf * math.cos(thetas[1])) * self.sin30
        x_2 = y_2 * self.tan60
        z_2 = -self.rf * math.sin(thetas[1])

        y_3 = (t + self.rf * math.cos(thetas[2])) * self.sin30
        x_3 = -y_3 * self.tan60
        z_3 = -self.rf * math.sin(thetas[2])

        dnm = (y_2 - y_1) * x_3 - (y_3 - y_1) * x_2

        w_1 = y_1 ** 2 + z_1 ** 2
        w_2 = x_2 ** 2 + y_2 ** 2 + z_2 ** 2
        w_3 = x_3 ** 2 + y_3 ** 2 + z_3 ** 2

        # x = (a1*z + b1)/dnm
        a_1 = (z_2 - z_1) * (y_3 - y_1) - (z_3 - z_1) * (y_2 - y_1)
        b_1 = -((w_2 - w_1) * (y_3 - y_1) - (w_3 - w_1) * (y_2 - y_1)) / 2.0

        # y = (a2*z + b2)/dnm
        a_2 = -(z_2 - z_1) * x_3 + (z_3 - z_1) * x_2
        b_2 = ((w_2 - w_1) * x_3 - (w_3 - w_1) * x_2) / 2.0

        # a*z^2 + b*z + c = 0
        a = a_1 ** 2 + a_2 ** 2 + dnm ** 2
        b = 2 * (a_1 * b_1 + a_2 * (b_2 - y_1 * dnm) - z_1 * dnm ** 2)
        c = (b_2 - y_1 * dnm) * (b_2 - y_1 * dnm) + b_1 * b_1 + dnm ** 2 * (z_1 ** 2 - self.re ** 2)

        # дискриминант
        d = b ** 2 - 4.0 * a * c
        if d < 0:
            return 0  # несуществующая позиция

        self.z = -0.5 * (b + math.sqrt(d)) / a
        self.x = (a_1 * self.z + b_1) / dnm
        self.y = (a_2 * self.z + b_2) / dnm
        return 1

    # обратная кинематика
    # вспомогательная функция, расчет угла theta1 (в плоскости YZ)
    def calculate_angle_yz(self, x_0, y_0, z_0):
        y_1 = -0.5 * 0.57735 * self.f  # f/2 * tg 30
        y_0 -= 0.5 * 0.57735 * self.e  # сдвигаем центр к краю
        
        # z = a + b*y
        a = (x_0 ** 2 + y_0 ** 2 + z_0 ** 2 + self.rf ** 2 - self.re ** 2 - y_1 ** 2) / (2 * z_0)
        b = (y_1 - y_0) / z_0

        # дискриминант
        d = -(a + b * y_1) * (a + b * y_1) + self.rf * (b ** 2 * self.rf + self.rf)
        if d < 0:
            return 0, None  # несуществующая точка
    
        y_j = (y_1 - a * b - math.sqrt(d)) / (b ** 2 + 1)  # выбираем внешнюю точку
        z_j = a + b * y_j
    
        theta = math.atan(-z_j / (y_1 - y_j)) * 180.0 / math.pi
        if y_j > y_1:
            theta += 180
        return 1, theta

    def inverse_kinematics(self, x_0, y_0, z_0):
        """
        обратная кинематика: (x0, y0, z0) -> (theta1, theta2, theta3)
        возвращаемый статус: 0=OK, -1=несуществующая позиция
        """
        status, self.theta_1 = self.calculate_angle_yz(x_0, y_0, z_0)  # first angle
        if status == 1:
            x_1 = x_0 * self.cos120 + y_0 * self.sin120
            y_1 = y_0 * self.cos120 - x_0 * self.sin120
            z_1 = z_0
            status, self.theta_2 = self.calculate_angle_yz(x_1, y_1, z_1)  # rotate coordinates to +120 deg
        if status == 1:
            x_2 = x_0 * self.cos120 - y_0 * self.sin120
            y_2 = y_0 * self.cos120 + x_0 * self.sin120
            z_2 = z_0
            status, self.theta_3 = self.calculate_angle_yz(x_2, y_2, z_2)  # rotate coordinates to -120 deg
        return status

    def print_angles(self):
        print("theta_1 = {}".format(self.theta_1))
        print("theta_2 = {}".format(self.theta_2))
        print("theta_3 = {}".format(self.theta_3))

    def print_coordinates(self):
        print("x = {}".format(self.x))
        print("y = {}".format(self.y))
        print("z = {}".format(self.z))


# размеры робота
# (обозначения см. на схеме)
len_big_base = 115.0
len_small_base = 457.3
len_big_shoulder = 232.0
len_small_shoulder = 112.0

TripodSolver = DeltaKinematics(len_big_base, len_small_base, len_big_shoulder, len_small_shoulder)
TripodSolver.forward_kinematics([100, 45, 45])
if not TripodSolver.inverse_kinematics(TripodSolver.x, TripodSolver.y, TripodSolver.z):
    print("Заданная позиция недостижима.")
else:
    TripodSolver.print_coordinates()
    TripodSolver.print_angles()
