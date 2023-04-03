from typing import Tuple
import numpy as np
from periodic_box import Box


def mass_center(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> Tuple[float, float, float]:
    if not (len(x) == len(y) == len(z)):
        raise ValueError(" Different length of arrays")
    x_m = np.mean(x)
    y_m = np.mean(y)
    z_m = np.mean(z)
    return x_m, y_m, z_m


def gyration_tensor(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    if not (len(x) == len(y) == len(z)):
        raise ValueError(" Different length of arrays")
    x_m, y_m, z_m = mass_center(x, y, z)
    x -= x_m
    y -= y_m
    z -= z_m
    T = np.zeros(shape=(3, 3), dtype=float)
    for i in range(len(x)):
        T[0][0] += x[i]*x[i]
        T[0][1] += x[i]*y[i]
        T[0][2] += x[i]*z[i]
        T[1][1] += y[i]*y[i]
        T[1][2] += y[i]*z[i]
        T[2][2] += z[i]*z[i]
    T[1][0] = T[0][1]
    T[2][0] = T[0][2]
    T[2][1] = T[1][2]
    return T/len(x)


def full_gyration_tensor(x: np.ndarray, y: np.ndarray, z: np.ndarray, box: Box) -> np.ndarray:
    if not (len(x) == len(y) == len(z)):
        raise ValueError(" Different length of arrays")
    T = np.zeros(shape=(3, 3), dtype=float)
    for i in range(len(x)):
        for j in range(len(x)):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dz = z[i] - z[j]
            dx, dy, dz = box.periodic_correct(dx, dy, dz)
            T[0][0] += dx*dx
            T[0][1] += dx*dy
            T[0][2] += dx*dz
            T[1][1] += dy*dy
            T[1][2] += dy*dz
            T[2][2] += dz*dz
    T[1][0] = T[0][1]
    T[2][0] = T[0][2]
    T[2][1] = T[1][2]
    return T/(2*len(x)**2)




class Tensor():
    def __init__(self, T):
        self.T = T
        self.Rx, self.Ry, self.Rz = sorted(np.linalg.eigvals(self.T))

    def __str__(self):
        return f"Rx = {self.Rx}, Ry = {self.Ry}, Rz = {self.Rz}"

    @property
    def principal_moments(self):
        return self.Rx, self.Ry, self.Rz

    @property
    def gyration_radius(self):
        return np.sqrt(self.Rx + self.Ry + self.Rz)

    @property
    def momets_ratio(self):
        return self.Rx / self.Rz, self.Ry / self.Rz, self.Ry / self.Rx

    @property
    def asphericity(self):
        return self.Rz - 0.5 * (self.Rx + self.Ry)

    @property
    def acylindricity(self):
        return self.Ry - self.Rx

    @property
    def shape_anisotropy(self):
        return 1.5 * (self.Rx**2 + self.Ry**2 + self.Rz**2)/((self.Rx + self.Ry + self.Rz)**2) - 0.5


if __name__ == "__main__":
    np.random.seed(42)
    # box_size = 60
    # box = Box(box_size, box_size, box_size)
    # r = 3  # radius of sphere
    # length = 0
    # center = (-3, -3, -3)
    # x, y, z = [], [], []
    # while length < 500:
    #     x_t, y_t, z_t = np.random.uniform(-r, r, 3)
    #     y_t = 0.
    #     z_t = 0.
    #     if (x_t**2 + y_t**2 + z_t**2 < r**2):
    #         x_t = x_t + center[0]
    #         y_t = y_t + center[1]
    #         z_t = z_t + center[2]
    #         # x_t, y_t, z_t = box.periodic_correct(x_t, y_t, z_t)
    #         length += 1
    #         x.append(x_t)
    #         y.append(y_t)
    #         z.append(z_t)
    # x = np.array(x)
    # y = np.array(y)
    # z = np.array(z)
    # m_center = mass_center(x, y, z)
    # g_tensor = gyration_tensor(x, y, z)
    # full_g_tensor = full_gyration_tensor(x, y, z, box)
    # # print(f"m_center = {m_center}")
    # # print(f"g_tensor = {g_tensor}")
    # # print(f"python_metho{np.linalg.eigvals(g_tensor)}")
    # # print(gyration_radius(x, y, z))
    # # print(np.sum(np.linalg.eigvals(g_tensor)))
    # # print(np.sum(np.linalg.eigvals(full_g_tensor)))
    # tensor = Tensor(g_tensor)
    # print(tensor)
    # print(tensor.acylindricity)
    # print(tensor.asphericity)
    # print(tensor.gyration_radius)
    # print(tensor.shape_anisotropy)
    # # plt.scatter(x, y)
    # # plt.show()
