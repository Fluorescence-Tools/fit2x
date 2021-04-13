@property
def x(self):
    return self.get_x()

@x.setter
def x(self, v):
    self.set_x(v)

@property
def y(self):
    return self.get_y()

@y.setter
def y(self, v):
    self.set_y(v)

@property
def ey(self):
    return self.get_ey()

@property
def squared_weights(self):
    return self.get_squared_weights()


