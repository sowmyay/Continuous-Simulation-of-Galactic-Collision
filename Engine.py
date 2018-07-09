
import numpy as np

class MainGalaxy(object):
    def __init__(self, distance_btw_galaxies, mass, z_gap=0, dense=False):
        self.distance = distance_btw_galaxies
        self.mass = mass
        self.stars = self.starPositions(distance_btw_galaxies, z_gap=z_gap, dense=dense)

    #Generates rings of stars for the main galaxy. Returns the total list of stars, specified by their position and velocity along x, y, z directions
    def starPositions(self, R, z_gap= 0 , dense=False):
        # star_count = 120
        stars_in_orbit = [12, 18, 24, 30, 36]
        layers = [5, 4, 3, 2, 1]
        radii = [20, 30, 40, 50, 60] if not dense else [12, 18, 24, 30, 36]
        stars = []

        # stars = pd.Dataframe(columns=["Mass", "PositionX", "PositionY", "PositionZ", "VelocityX", "VelocityY", "VelocityZ"])
        for i, r in enumerate(radii):
            count = stars_in_orbit[i]
            theta_diff = 360 / count
            theta = 0
            z_gap = z_gap
            G = 4.302 * np.exp(-3)
            w = np.sqrt(G * self.mass / r)
            for index, j in enumerate(range(layers[i])):
                s = self.one_ring(z=z_gap * index, R=R, r=r, w=w, theta=theta, theta_diff=theta_diff, stars_in_orbit=stars_in_orbit[i])
                stars += s

        return stars

    # Generates a single ring of stars, and returns an array of the position of each of the star.
    # It takes in the z position as input, and constructs a ring of stars in a plane parallel to the xy plane
    def one_ring(self, z, R, r, w, theta, theta_diff, stars_in_orbit):
        stars = []
        for j in range(stars_in_orbit):
            radius = R * r / 100
            star = [1, radius * np.cos(theta), radius * np.sin(theta), z, radius * w * np.cos(theta - 90),
                    radius * w * np.sin(theta - 90), 0]
            stars += [star]
            theta += theta_diff
        return stars

class DisruptorGalaxy(object):
    def __init__(self, M, initial_conditions):
        self.mass = M
        # self.velocity = v
        # self.theta =theta*0.0174533
        # self.phi = phi*0.0174533
        self.R = initial_conditions
        # self.R = [distance*np.cos(theta)*np.cos(phi), distance*np.sin(theta)*np.cos(phi), distance*np.sin(phi), v*np.cos(theta)*np.cos(phi), v*np.sin(theta)*np.cos(phi), v*np.sin(phi)]
        # self.R = [-8, 9, -9, 0.85, -0.5, -0.2]

class Engine(object):
    def __init__(self,  M = 330, S=330 , initial_conditions=[-8, -9, 0, 0.85, 0.65, 0], z_gap=2, dense=False): #distance in parsecs
        posx = initial_conditions[0]
        posy = initial_conditions[1]
        posz = initial_conditions[2]
        distance = np.sqrt(pow(posx, 2)+ pow(posy, 2)+ pow(posz, 2))
        self.main = MainGalaxy(distance, M, z_gap=z_gap, dense=dense) #rad/s
        self.disruptor = DisruptorGalaxy(S, initial_conditions=initial_conditions)









