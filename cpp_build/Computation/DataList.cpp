#include <vector>
#include <iostream>

#include "State.h"
#include "DState.h"
#include "DataList.h"
#include "pairwise_connections.h"



// Populate datalist with states and mesh connectivity data
template <typename T>
DataList<T>::DataList()
{

}

template <typename T>
DataList<T>::DataList(std::vector<T> data, std::vector<Pairwise_Connection> connectivity, std::vector<Pairwise_Rigid_Connection> rig_connectivity)
{
    this->data = data;
    this->connectivity = connectivity;
    this->rig_connectivity = rig_connectivity;
}

//make a copy of a given datalist (not just reference it in memory) and return the copy
template <typename T>
DataList<T> DataList<T>::clone()
{
    return  DataList(this->data,this->connectivity,this->rig_connectivity);
}

template <typename T>
DataList<T> DataList<T>::combine(DataList<T> datalist)
{
    std::vector<T> tempsarr;
    std::vector<Pairwise_Connection> tempcarr;
    std::vector<Pairwise_Rigid_Connection> temprarr;

    // Combine state data
    for (int i = 0; i < this->data.size(); i++)
    {
        tempsarr.push_back(this->data[i]);
    }
    for (int i = 0; i < datalist.data.size(); i++)
    {
        tempsarr.push_back(datalist.data[i]);
    }

    // Combine connectivity data
    if (this->connectivity.size() != 0)
    {
        for (int i = 0; i < this->connectivity.size(); i++)
        {
            tempcarr.push_back(this->connectivity[i]);
        }
    }
    if (datalist.connectivity.size() != 0)
    {
        for (int i = 0; i < datalist.connectivity.size(); i++)
        {
            Pairwise_Connection temp_pair;
            // Base vertices for relabeling when adding new set of vertices
            int vert_num = this->data.size();
            temp_pair.index_pair[0] = datalist.connectivity[i].index_pair[0] + vert_num;
            temp_pair.index_pair[1] = datalist.connectivity[i].index_pair[1] + vert_num;
            tempcarr.push_back(temp_pair);
        }
    }

    // Combine rig_connectivity data
    if (this->rig_connectivity.size() != 0)
    {
        for (int i = 0; i < this->rig_connectivity.size(); i++)
        {
            temprarr.push_back(this->rig_connectivity[i]);
        }
    }
    if (datalist.rig_connectivity.size() != 0)
    {
        for (int i = 0; i < datalist.rig_connectivity.size(); i++)
        {
            Pairwise_Rigid_Connection temp_pair;
            // Base vertices for relabeling when adding new set of vertices
            int vert_num = this->data.size();
            temp_pair.index_pair[0] = datalist.rig_connectivity[i].index_pair[0] + vert_num; 
            temp_pair.index_pair[1] = datalist.rig_connectivity[i].index_pair[1] + vert_num;
            temp_pair.lagrange_multipler = datalist.rig_connectivity[i].lagrange_multipler;
            temprarr.push_back(temp_pair);
        }
    }
    return DataList<T>(tempsarr,tempcarr,temprarr);
    
}

template <>
std::vector<double> DataList<State>::toArray()
{
    std::vector<double> temparr;

    for (int i = 0; i < this->data.size(); i++)
    {
        temparr.insert(std::end(temparr), std::begin(this->data[i].pos), std::end(this->data[i].pos));
        temparr.insert(std::end(temparr), std::begin(this->data[i].vel), std::end(this->data[i].vel));
    }

    return temparr;

}

template <>
std::vector<double> DataList<DState>::toArray()
{
    std::vector<double> temparr;

    for (int i = 0; i < this->data.size(); i++)
    {
        temparr.insert(std::end(temparr), std::begin(this->data[i].vel), std::end(this->data[i].vel));
        temparr.insert(std::end(temparr), std::begin(this->data[i].acc), std::end(this->data[i].acc));
    }

    return temparr;

}

//print information about the state to the terminal
template<>
void DataList<State>::print_info()
{
    int count = 0;
    // Add the elements of state.vec and this->vel
    for (State item : this->data)  
    {  
        std::cout << "State " << count << " :\n";
        item.print_info();  
        count++;
    };
}

template<>
void DataList<DState>::print_info()
{
    int count = 0;
    // Add the elements of state.vec and this->vel
    for (DState item : this->data)  
    {  
        std::cout << "DState " << count << " :\n";
        item.print_info();  
        count++;
    };
}

template <typename T>
void DataList<T>::add(DataList<T> datalist)
{
    // Add the elements of state.vec and this->vel
    for (int i = 0; i < this->data.size(); i++)  
    {  
        T temp = datalist.data[i];
        this->data[i].add(temp);  
    };
}

template <typename T>
void DataList<T>::sub(DataList<T> datalist)
{
    // Add the elements of state.vec and this->vel
    for (int i = 0; i < this->data.size(); i++)  
    {  
        T temp = datalist.data[i];
        this->data[i].sub(temp);  
    };
}

template <typename T>
void DataList<T>::multiplyScalar(double k)
{
    // Add the elements of state.vec and this->vel
    for (int i = 0; i < this->data.size(); i++)  
    {  
        this->data[i].multiplyScalar(k);  
    };
}

template <>
void DataList<State>::flow(double eps)
{
    // Add the elements of state.vec and this->vel
    for (State item : this->data)  
    {  
        item.flow(eps);  
    };
}

    // # //implementing .updateBy() componentwise
    // # //WARNING!!! RIGHT NOW NOT ALL OBJECTS HAVE AN UPDATEBY
    // # //ONLY STATE, AND IT TAKES DSTATE: SO THIS WILL RETURN
    // # //AN ERROR IF ITS CALLED IN A CASE WHERE THE INDIVIDUAL OBJECTS DONT
    // # //IMPLEMENT IT:  THAT'S EXPECTED; BUT WE SHOULD INCLUDE AN ACTUAL ERROR
    // # //MESSAGE TO CONSOLE.LOG() THAT EXPLAINS.


template <>
void DataList<State>::updateBy(DataList<DState> datalist)
{
// Add the elements of state.vec and this->vel
for (int i = 0; i < this->data.size(); i++)  
{  
    DState temp = datalist.data[i];
    this->data[i].updateBy(temp);  
};
}

// //The explicit instantiation part
template class DataList<State>; 
template class DataList<DState>;