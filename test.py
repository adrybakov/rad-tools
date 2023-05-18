class A:
    def __init__(self, a, b) -> None:
        self.a = a
        self.b = b

    @property
    def c(self):
        return self.a + self.b


class B(A):
    def __init__(self, a, b, c) -> None:
        self.a = a
        self.b = b
        self.c = c


a = A(1, 2)
b = B(1, 2, 4)
print(a.a, a.b, a.c)
print(b.a, b.b, b.c)
