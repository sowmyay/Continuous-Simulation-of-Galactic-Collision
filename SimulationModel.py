
import numpy as np
import pandas as pd
from Engine import Engine
from matplotlib import pyplot
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D


G = 4.302 * np.exp(-3) #Gravitational constant for when
                        # distances are measured in parsecs (pc),
                        # velocities in kilometers per second (km/s)
                        # and masses in solar units


class SimulationApplication(object):


    def __init__(self, masses, engine_inputs, time = 20, time_steps=12):
        # time_steps = pow(10, 12)

        self.engine = Engine(M =masses[0], S=masses[1] , initial_conditions=engine_inputs[:6], z_gap=engine_inputs[6], dense=engine_inputs[7])
        new_stars, dis = self.runSim(time=time, timesteps=time_steps)
        self.plotNewStars(new_stars, dis,
                          stars_count=len(self.engine.main.stars),
                          timesteps=time_steps)

    # computes updated position for the disruptor galaxy and every star of the main galaxy at the specified number of time steps, upto the total time
    # returns list of new_starsx, new_stary, new_starsz, which are dataframes of new x, y, z positions at every time step. - rows = star, columns = range(time_steps)
    #               and  updated postion dataframe of the disruptor galaxy and stars of main galaxy in a dataframe - rows = (x, y, z), columns = range(time_steps)

    def runSim(self, time=20, timesteps=2):
        stars = self.engine.main.stars
        disrup = self.engine.disruptor
        M = self.engine.main.mass
        S = self.engine.disruptor.mass

        t0 = 0
        dis = disrup.R
        new_starsx = pd.DataFrame(np.nan, index=range(len(stars)), columns=range(timesteps)) #empty dataframe to fill the updated x postion of each star(row) at each timestep (column)
        new_starsy = pd.DataFrame(np.nan, index=range(len(stars)), columns=range(timesteps)) #empty dataframe to fill the updated y postion of each star(row) at each timestep (column)
        new_starsz = pd.DataFrame(np.nan, index=range(len(stars)), columns=range(timesteps)) #empty dataframe to fill the updated z postion of each star(row) at each timestep (column)
        disruptor = pd.DataFrame(np.nan, index=range(3), columns=range(timesteps)) #empty dataframe to fill the updated x, y, z postion of disruptor galaxy(rows) at each timestep (column)
        locs = stars[:]
        jump = np.log10(time)/timesteps
        for j in range(timesteps):
            t1 = pow(10, (j+1)*jump) # time at the end of each time step jump
            dis = self.evalstar(func=self.disrode,
                                y0=np.array(dis), M=M, S=S, t0=t0, t1=t1,
                                Rv=disrup.R,
                                timesteps=2).tolist()[1]
            disruptor.iloc[0, j] = dis[0]
            disruptor.iloc[1, j] = dis[1]
            disruptor.iloc[2, j] = dis[2]

            for i, s in enumerate(stars):
                print i
                y = self.evalstar(func=self.starode,
                                  y0=np.array(locs[i]), M=M, S=S, t0=t0, t1=t1,
                                  Rv=dis,
                                  timesteps=2)
                locs[i] = y.tolist()[1]
                temp = pd.DataFrame(y)
                temp = temp.iloc[:, 1:4]
                temp = temp.transpose()
                new_starsx.iloc[i, j] = temp.iloc[0, 0]
                new_starsy.iloc[i, j] = temp.iloc[1, 0]
                new_starsz.iloc[i, j] = temp.iloc[2, 0]
        return [new_starsx, new_starsy, new_starsz], disruptor

    #Plots the positions of the stars of the main galaxy and disruptor galaxy at the specified number of time steps, and saves it as image files
    def plotNewStars(self, new_stars, disruptor, stars_count, timesteps):
        new_starsx = new_stars[0]
        new_starsy = new_stars[1]
        new_starsz = new_stars[2]
        for t in range(timesteps):
            fig = pyplot.figure()
            ax = Axes3D(fig)
            x_pos = disruptor.iloc[0, t]
            y_pos = disruptor.iloc[1, t]
            z_pos = disruptor.iloc[2, t]
            ax.scatter(x_pos, y_pos, z_pos, c="blue")
            for i in range(stars_count):
                x_pos = new_starsx.iloc[i, t]
                y_pos = new_starsy.iloc[i, t]
                z_pos = new_starsz.iloc[i, t]
                ax.scatter(x_pos, y_pos, z_pos, c="red")

            fig.savefig("Timestep_" + str(t))
            # pyplot.show()

    # computes updated position for a single star of the main galaxy
    # returns dy = [0, vx, vy, vz, ax, ay, az]
    #               |____of disruptor_______|
    def starode(self,  y, t, M, S, Rv):
        r = max(np.sqrt(y[1] ** 2 + y[2] ** 2 + y[3] ** 2), 0.1)
        # new_R = Rv
        # Rv[0] = Rv[0] - Rv[3] * t # * np.cos(theta) * np.cos(phi)
        # Rv[1] = Rv[1] - Rv[4] * t #* np.sin(theta) * np.cos(phi)
        # Rv[2] = Rv[2] - Rv[5] * t #* np.sin(phi)

        R = max(np.sqrt(Rv[0] ** 2 + Rv[1] ** 2 + Rv[2] ** 2), 0.1)
        pv = np.subtract(Rv[:3], y[1:4])
        p = np.sqrt(pv[0] ** 2 + pv[1] ** 2 + pv[2] ** 2)
        dy = [0]
        dy += [y[4]]  # vx
        dy += [y[5]]  # vy
        dy += [y[6]]  # vz
        dy += [-1*G * ((M / pow(r, 3)) * y[1] - (S / pow(p, 3)) * pv[0] + (S / pow(R, 3)) * Rv[0])]  # ax
        dy += [-1*G * ((M / pow(r, 3)) * y[2] - (S / pow(p, 3)) * pv[1] + (S / pow(R, 3)) * Rv[1])]  # ay
        dy += [-1*G * ((M / pow(r, 3)) * y[3] - (S / pow(p, 3)) * pv[2] + (S / pow(R, 3)) * Rv[2])] # az
        return dy

    # computes updated position for the disruptor galaxy
    # returns dy = [vx, vy, vz, ax, ay, az]
    #               |____of disruptor___|
    def disrode(self, y, t, M, S, Rv):

        Rv[0] = Rv[0] - Rv[3] * t # * np.cos(theta) * np.cos(phi)
        Rv[1] = Rv[1] - Rv[4] * t #* np.sin(theta) * np.cos(phi)
        Rv[2] = Rv[2] - Rv[5] * t #* np.sin(phi)

        R = max(np.sqrt(Rv[0]**2+Rv[1]**2+Rv[2]**2), 0.1)
        # dy = [0]
        dy = []
        dy += [y[3]]  # vx
        dy += [y[4]]  # vy
        dy += [y[5]]  # vz
        dy += [-1*G*(M+S)/R**3*Rv[0]] #ax
        dy += [-1*G*(M+S)/R**3*Rv[1]] #ay
        dy += [-1*G*(M+S)/R**3*Rv[2]] #az
        return dy

    # Calls the ode function specific to the disruptor or the star,
    # and then uses SciPy function odeint() to solve the system of odes using Runge-Kutta numerical methods for integration
    # returns y = [x, y, z, vx, vy, z] or [0, x, y, z, vx, vy, z]
    #              |__if disruptor___|    |__if star___________|
    def evalstar(self, func, y0, M, S, Rv, t0=0, t1=pow(10, 3), timesteps=2):
        t = np.linspace(t0, t1, timesteps)  # points of evaluation of solution
        y = np.zeros((len(t), len(y0)))  # array for solution
        y[0, :] = y0

        for i in range(1, t.size):
            # print y0
            y = odeint(func, t=t, y0=y0, args=(M, S, Rv)) ## Solve the differential equation for each time index
                                         ##and output the position values of the stars and disruptor galaxy.
        return y


if __name__ == "__main__":
    inp = np.loadtxt("input.txt")
    initial_conditions = [-8, -9, 0, 0.85, 0.65, 0, 2, 0]
    # -8  # posX
    # -9  # posY
    # 0  # posZ
    # 0.85  # vx
    # 0.65  # vy
    # 0  # vz
    # 2  # z_gap
    # 0  # dense
    s = SimulationApplication(masses=inp[:2], engine_inputs=initial_conditions, time = inp[2], time_steps=int(inp[3]))
