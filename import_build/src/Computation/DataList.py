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
        array of state objects

    Methods
    -------
    clone()
        Generate copy of self
        Returns DataList clone

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
    def __init__(self, data):
        # data is list of State objects (describing motion of vertices)
        self.data = data

        # //implementation of .clone() for the array
    def clone(self):
        temparr = []
        for a in self.data:
            temparr.append(a.clone())
        return  DataList(temparr)

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