class grid:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.grid = [[0 for x in range(width)] for y in range(height)]
    def __str__(self):
        return str(self.grid)
    def replace(self, x, y, value):
        if value > self.grid[y][x]:
            self.grid[y][x] = value
    def set(self, x, y, value):
        self.grid[y][x] = value
    def get(self, x, y):
        return self.grid[y][x]
    def get_width(self):
        return self.width
    def get_height(self):
        return self.height
    def get_grid(self):
        return self.grid
    def get_neighbours(self, x, y):
        neighbours = []
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i == 0 and j == 0:
                    continue
                if x + i < 0 or x + i >= self.width:
                    continue
                if y + j < 0 or y + j >= self.height:
                    continue
                neighbours.append((x + i, y + j))
        return neighbours
    def get_neighbours_value(self, x, y):
        neighbours = []
        for i in range(-1, 2):
            for j in range(-1, 2):
                if i == 0 and j == 0:
                    continue
                if x + i < 0 or x + i >= self.width:
                    continue
                if y + j < 0 or y + j >= self.height:
                    continue
                neighbours.append(self.grid[y + j][x + i])
        return neighbours

if __name__ == '__main__':
    g = grid(10, 10)
    g.set(5, 5, 100)
    g.replace(5, 5, 60)
    print(g)

