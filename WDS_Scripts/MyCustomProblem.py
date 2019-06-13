from platypus import NSGAII, Problem, Integer,Real
import sys
sys.path.insert(0, '/home/varsha/Documents/MyCode/Water Network/Optimization/')
import Functions
import Settings
class my_mo_problem(Problem):
	et, hStar, o_curves, n_curves, nbOfPipes,nbOfPumps,nbOfTanks,Conn,NoConn,max_elevation = Settings.SetValues('output/WCR-modified-12018_11_01__22_05_55.inp')
	Functions.SetVariables(et)
	def __init__(self):
		super(my_mo_problem, self).__init__(self.nbOfPipes + 25*self.nbOfPumps + 3*self.nbOfTanks ,2, 1)
		self.types[:] = [Real(0, 9)]*self.nbOfPipes + [Real(0,self.n_curves-1)]*self.nbOfPumps + [Real(25,100)]*self.nbOfTanks + [Real(25,40)]*self.nbOfTanks+ [Real(9,10)]*self.nbOfTanks+[Real(0,1)]*(24*self.nbOfPumps)
		self.constraints[:] = "<=0"
		self.directions[:] = Problem.MINIMIZE
	def evaluate(self, solution):
		pipes = solution.variables[0:self.nbOfPipes] #diameter of pipe
		pumps = solution.variables[self.nbOfPipes:self.nbOfPipes+self.nbOfPumps] #curve of pumps
		tanks_diam = solution.variables[self.nbOfPipes+self.nbOfPumps:self.nbOfPipes+self.nbOfPumps+self.nbOfTanks] #diameter of tank
		tanks_max =solution.variables[self.nbOfPipes+self.nbOfPumps+self.nbOfTanks:self.nbOfPipes+self.nbOfPumps+2*self.nbOfTanks] # max level of tank
		tanks_min = solution.variables[self.nbOfPipes+self.nbOfPumps+2*self.nbOfTanks:self.nbOfPipes+self.nbOfPumps+3*self.nbOfTanks] #min level of tank
		patterns = solution.variables[self.nbOfPipes + self.nbOfPumps + 3 * self.nbOfTanks:self.nbOfPipes + self.nbOfPumps + 3 * self.nbOfTanks + 24 * self.nbOfPumps]
		solution.objectives[:] = [-Functions.Res(pipes,patterns,pumps,tanks_diam,tanks_max,tanks_min,self.et,self.hStar,self.n_curves,self.Conn,self.NoConn,self.max_elevation),Functions.Cost(patterns,self.et)]
		solution.constraints[:] = Functions.Constraint()
