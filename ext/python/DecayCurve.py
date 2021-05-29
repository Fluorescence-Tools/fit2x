@property
def x(self):
    return np.array(self.get_x())

@x.setter
def x(self, v):
    self.set_x(v)

@property
def y(self):
    return np.array(self.get_y())

@y.setter
def y(self, v):
    self.set_y(v)

@property
def ey(self):
    return np.array(self.get_ey())

