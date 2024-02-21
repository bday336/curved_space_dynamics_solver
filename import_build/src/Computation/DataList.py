# // a class for working with points in the configuration space of a
# //composite system.
# //instances of the DataList class are arrays of objects
# //these objects need to implement the methods .clone(), .add(), .sub(), and .multiplyScalar()
# //thus, DataLists can be lists of either states, or dStates as required:
import numpy as np

class DataList:
    """
    A class used to store information about collection of states (collection of vertices)

    ...

    Attributes
    ----------
    data : array
        array of state objects (vertices)

    connectivity : array
        array of connectivity information between states in data - defaults to None i.e. data is collection of free vertices
        (e.g. connectivity = [[1,0],[3,8]] means states 0 and 1 are connected and states 3 and 8 are connected via coupling potential or constraint)

    Methods
    -------
    clone()
        Generate copy of self
        Returns DataList clone

    combine(dataList)
        Generates new DataList object consisting of states from self and dataList 
        Returns combined Datalist object

    add(dataList)
        Add velocity of each state in dataList to each state in self componentwise
        Returns self

    sub(dataList)
        Subtract velocity of each state in dataList from each state in self componentwise
        Returns self

    multiplyScalar(k)
        Scale each state in self by scalar value k
        Returns self

    flow(eps)
        Move each state in self infinitesimally (i.e. by eps) along its tangent vector
        Returns self

    updateBy(dataList)
        Update each state in self by infinitesimally flowing along its respective differential (dState) contained in dataList (DataList of dStates)
        Returns self
    """

    # Populate list with states
    def __init__(self, data, connectivity = []):
        # data is list of State objects (describing motion of vertices)
        self.data = data
        self.connectivity = connectivity

        # //implementation of .clone() for the array
    def clone(self):
        temparr = []
        for a in self.data:
            temparr.append(a.clone())
        return  DataList(temparr, self.connectivity)
    
    def combine(self, dataList):
        tempsarr = []
        tempcarr = []

        # Combine state data
        for a in self.data:
            tempsarr.append(a.clone())
        for b in dataList.data:
            tempsarr.append(b.clone())

        # Combine connectivity data
        if len(self.connectivity) != 0:
            tempcarr = tempcarr + self.connectivity
        if len(dataList.connectivity) != 0:
                for c in dataList.connectivity:
                    c = [c[0] + len(self.data), c[1] + len(self.data)]
                    tempcarr.append(c)
        return DataList(tempsarr,tempcarr)
    
    # Euclidean norm for implicit solver methods (Only applicable for DataLists of dStates currently)
    def norm(self):
        val = 0.
        for a in self.data:
            val = a.vel[0]**2. + a.vel[1]**2. + a.vel[2]**2. + a.acc[0]**2. + a.acc[1]**2. + a.acc[2]**2.

        val = np.sqrt(val)
        return val
    
    def toArray(self):
        res_array = []

        # For DataList of States
        if self.data[0].__class__.__name__ == "State":
            for c in self.data:
                res_array.append(c.pos.copy())
                res_array.append(c.vel.copy())

        # For DataList of dStates
        if self.data[0].__class__.__name__ == "dState":
            for c in self.data:
                res_array.append(c.vel.copy())
                res_array.append(c.acc.copy())

        res_array = np.array(res_array).flatten()

        return res_array
        

    # //implementing .add() componentwise (between datalists)
    def add(self, dataList ):
        for a in range(len(self.data)):
            # print(a)
            self.data[a] = self.data[a].add(dataList.data[a])
        return self

    # //implementing .sub() componentiwse (between datalists)
    def sub(self, dataList ):
        for a in range(len(self.data)):
            # print(a)
            self.data[a] = self.data[a].sub(dataList.data[a])
        return self
    
    # //implementing .multiplyScalar() componentwise
    def multiplyScalar(self, k ):
        for a in self.data:
            a = a.multiplyScalar(k)
        return self

    def flow(self, eps):
        for a in self.data:
            a = a.flow(eps)
        return self

    # //implementing .updateBy() componentwise
    # //WARNING!!! RIGHT NOW NOT ALL OBJECTS HAVE AN UPDATEBY
    # //ONLY STATE, AND IT TAKES DSTATE: SO THIS WILL RETURN
    # //AN ERROR IF ITS CALLED IN A CASE WHERE THE INDIVIDUAL OBJECTS DONT
    # //IMPLEMENT IT:  THAT'S EXPECTED; BUT WE SHOULD INCLUDE AN ACTUAL ERROR
    # //MESSAGE TO CONSOLE.LOG() THAT EXPLAINS.

    def updateBy(self, dataList ):
        for a in range(len(self.data)):
            # print(a)
            self.data[a] = self.data[a].clone().updateBy(dataList.data[a])
            # print(a)
            # print(self.data[a].pos)
        return self