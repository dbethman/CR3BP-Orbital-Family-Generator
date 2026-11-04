import asset_asrl as ast
from asset_asrl.Astro.AstroModels import CR3BP
import asset_asrl.Astro.Constants as c
import numpy as np



class FamilySolver:
    def __init__(self, mu1=c.MuEarth, mu2=c.MuMoon, lstar=c.LD, nSeg=800):
        self.oc = ast.OptimalControl
        self.vf = ast.VectorFunctions

        self.PhaseRegion = self.oc.PhaseRegionFlags
        self.TransModes = self.oc.TranscriptionModes

        self.mu1 = mu1
        self.mu2 = mu2
        self.lstar = lstar
        self.nSeg = nSeg

        self.dt = np.pi / 10000

        self.ode = CR3BP(self.mu1, self.mu2, self.lstar)

        self.mu = self.ode.mu

        self.odeIntegrator = self.ode.integrator(self.dt)

    def solve(self, init, tf, fixedvars=[0,1,2], dS=None, dV=None, V2=None):
        steps = 1000
        trajguess = self.odeIntegrator.integrate_dense(init, tf, steps)
        odePhase = self.ode.phase(self.TransModes.LGL3)

        

        odePhase.setTraj(trajguess, self.nSeg)

        if dS is None:
            # --- STANDARD SHOOTING CONSTRAINTS ---
            for idx in fixedvars:
                odePhase.addBoundaryValue(self.PhaseRegion.Front, [idx], [trajguess[0][idx]])

            odePhase.addBoundaryValue(
                    self.PhaseRegion.Front, [3, 5, 6], [0.0, 0.0, 0.0] 
                )
            odePhase.addBoundaryValue(
                    self.PhaseRegion.Back, [1, 3, 5], [0.0, 0.0, 0.0] 
                )
            
        else:
            # --- PAC CONSTRAINTS ---
            odePhase.addBoundaryValue(self.PhaseRegion.Front, [1, 3, 5, 6], [0.0, 0.0, 0.0, 0.0])
            odePhase.addBoundaryValue(self.PhaseRegion.Back, [1, 3, 5], [0.0, 0.0, 0.0])

            V2 = [float(v) for v in V2]
            dV = [float(v) for v in dV]

            args = self.vf.Arguments(14)
            X0 = args.head(7)
            XF = args.tail(7)

            pac = (dV[0] * (X0[0] - V2[0])
                 + dV[1] * (X0[2] - V2[1])
                 + dV[2] * (X0[3] - V2[2])
                 + dV[3] * (X0[4] - V2[3])
                 + dV[4] * (X0[5] - V2[4])
                 + dV[5] * (XF[6] - V2[5])
                 - dS)

            odePhase.addEqualCon(self.PhaseRegion.FrontandBack, pac, list(range(7)))
            odePhase.setThreads(10, 10)


        
        odePhase.optimizer.EContol = 1e-11
        odePhase.solve()

        lastiter = odePhase.optimizer.LastIterNum

        trajSol = odePhase.returnTraj()

        return trajSol, lastiter
    

    def continuation(self, init, tf, init_step_x,init_step_z, min_distance, dS=.001, dS_min=0.000001, maxiters=1000):

        trajectories = []
        
        
        trajSol, _ = self.solve(init, tf) 
        trajectories.append(trajSol)

        init = list(trajSol[0]) 
        tf = trajSol[-1][6]
        V1 = np.array([init[0], init[2], init[3], init[4], init[5], tf])
        init[0] += init_step_x
        init[2] += init_step_z

        trajSol, _ = self.solve(init, tf)
        trajectories.append(trajSol)

        init = list(trajSol[0])
        tf = trajSol[-1][6]
        V2 = np.array([init[0], init[2], init[3], init[4], init[5], tf])

        
        


        i = 0
        while i < maxiters:
            dV = V2 - V1
            dV = dV / np.linalg.norm(dV) 

            try:
                if dS < dS_min:
                    print(f"Continuation stopped: Minimum dS ({dS_min}) reached at step {i}.")
                    break
                trajSol, iters = self.solve(init, tf, dS=dS, dV=dV, V2=V2)

                if iters > 20:
                    dS *= 0.8
                    continue

                
                
                if iters < 5:
                    dS *= 1.2
                traj_arr = np.array(trajSol)
                positions = traj_arr[:, 0:3]

                earth_position  = np.array([-self.mu, 0.0, 0.0])
                moon_position = np.array([1.0 - self.mu, 0.0, 0.0])

                earth_distance = np.linalg.norm(positions - earth_position , axis=1)
                moon_distance = np.linalg.norm(positions - moon_position, axis=1)
                closest_approach = min(np.min(earth_distance), np.min(moon_distance))

                

                
                trajectories.append(trajSol) 

                init = list(trajSol[0]) 
                tf = trajSol[-1][6]
                V1 = V2
                V2 = np.array([init[0], init[2], init[3], init[4], init[5], tf])


                arr = np.array(trajSol)

                if closest_approach < min_distance:
                                    print(f"Continuation stopped: Closest approach ({closest_approach:.5f}) below minimum distance ({min_distance:.5f}) at step {i}.")
                                    break
                i+=1
 
            except Exception as e:
                dS *= 0.7
                if dS < dS_min:
                    print(f"Continuation stopped: Solver failed and minimum dS ({dS_min}) reached at step {i}.")
                    break
                continue




            
        return trajectories



