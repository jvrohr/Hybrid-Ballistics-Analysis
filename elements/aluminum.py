class Aluminum:
    A = [4.8, 0.00322, 155.239]

    def get_CP(self, temperature: float):
        return (self.A[0] + self.A[1]*temperature)*self.A[2]