#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }

    //calculate the _maxPinNum
    for(int i = 0; i < _cellArray.size(); ++i)
    {
        if(_cellArray[i]->getPinNum() > _maxPinNum)
        {
            _maxPinNum = _cellArray[i]->getPinNum();
        }
    }
    //initialize the bucket list
    _bl[0] = new Bucketlist(_maxPinNum);
    _bl[1] = new Bucketlist(_maxPinNum);

    //set up the lower and upper bound of the cell number in each partition
    _min_cell_num = floor((1 - _bFactor) * _cellNum / 2) + 1;
    _max_cell_num = ceil((1 + _bFactor) * _cellNum / 2) - 1;

    // printf("min_cell_num = %d\n", _min_cell_num);
    // printf("max_cell_num = %d\n", _max_cell_num);


    return;
}

void Partitioner::partition()
{
    initPartition();
    cutSize();
    
    //initialize the params
    _maxAccGain = 0;
    _accGain = 0;
    _iterNum = 0;
    _moveNum = 0;
    _bestMoveNum = 0;
    
    //test below
    // printSummary();
    // _moveStack.clear();
    // initCellGain();
    // checkCellNodePointers();
    // unlockAllCell();
    // int max_gain_cell_id = findMaxGainCell();
    // Cell* max_gain_cell = _cellArray[max_gain_cell_id]; 
    // _moveStack.push_back(max_gain_cell_id);

    // _accGain += max_gain_cell->getGain();
    // _maxAccGain = max(_maxAccGain, _accGain);
    // _bestMoveNum = max(_bestMoveNum, _moveNum);

    // max_gain_cell->lock();
    // updateCellGain(max_gain_cell_id);
    // moveCell(max_gain_cell_id);
    // cutSize();
    // printSummary();
    //test above

    do
    {
        _iterNum++;
        //printf("================iterNum = %d====================\n", _iterNum);
        //initialize the cell gain
        _moveStack.clear();
        unlockAllCell();

        initCellGain();
        // printSummary();


        _maxAccGain = 0;
        _accGain = 0;
        _bestMoveNum = 0;
        _moveNum = 0;
        
        do
        {
            _moveNum++;
            // printf("===============moveNum = %d=================\n", _moveNum);
            // printf("partSize[0] = %d\n", _partSize[0]);
            // printf("partSize[1] = %d\n", _partSize[1]);
            int max_gain_cell_id = findMaxGainCell();

            //check the max gain cell id
            if(max_gain_cell_id == -1)
            {
                printf("error: max gain cell id is -1\n");
                continue;
            }
            // printf("max_gain_cell_id = %d\n", max_gain_cell_id);
            bool from = _cellArray[max_gain_cell_id]->getPart();
            bool to = !from;
            // // //check if the cell is locked
            // printf("_bl[0] bucketlist\n");
            // _bl[0]->printBucketlist();
            // printf("\n_bl[1] bucketlist\n");
            // _bl[1]->printBucketlist();


            //delete the max gain node from the bucket list
            // cout << "deletion of max gain node" << endl;
            _bl[from]->deleteNode(_cellArray[max_gain_cell_id]->getNode());
            // cout << "end of deletion of max gain node" << endl;
            //illegal move because the two partitions are out of bound

            Cell* max_gain_cell = _cellArray[max_gain_cell_id];
            
            _moveStack.push_back(max_gain_cell_id);


            //params update
            _accGain += max_gain_cell->getGain();
            // printf("accGain = %d\n", _accGain);
            if(_accGain > _maxAccGain)
            {
                _maxAccGain = _accGain;
                _bestMoveNum = _moveNum;
            }


            //FM heuristic
            max_gain_cell->lock();
          

            // cout << "update gain for cell " << endl;
            updateCellGain(max_gain_cell_id);
            // cout << "end of update gain for cell " << endl;

            moveCell(max_gain_cell_id);
            // cutSize();
            // printf("cutSize = %d\n", _cutSize);

            
            
        }while(_moveNum < _cellNum);
        //printf("maxAccGain = %d\n", _maxAccGain);
        //printf("bestMoveNum = %d\n", _bestMoveNum);
        retreatOptimal();
        //cutSize();
        //printf("retreat optimal, result cutSize = %d\n", _cutSize);


    } while (_maxAccGain > 0);
    
    


 
}


void Partitioner::retreatOptimal()
{
    for(int i = _bestMoveNum; i < _moveStack.size(); ++i)
    {
        moveCell(_moveStack[i]);
    }
}

void Partitioner::checkCellNodePointers() const
{
    printf("=== Debug: Checking each cell's node pointers ===\n");
    for (size_t i = 0; i < _cellArray.size(); ++i)
    {
        Node* node = _cellArray[i]->getNode();
        if (node == nullptr)
        {
            printf(" Cell ID %zu has null node pointer!\n", i);
            continue;
        }

        if (node->getPrev() == nullptr)
        {
            printf(" Cell ID %d | Node ID %d has NULL prev pointer!\n", (int)i, node->getId());
        }

        // Optional: Show full node info
        printf("Cell ID %d | Node ID %d | Prev: %p | This: %p | Next: %p\n",
               (int)i,
               node->getId(),
               (void*)node->getPrev(),
               (void*)node,
               (void*)node->getNext());
    }
    printf("=== End of Cell Node Pointer Check ===\n");
}




int Partitioner::findMaxGainCell()
{
    // bool legal_0 = (_partSize[0] > _min_cell_num && _partSize[0] < _max_cell_num);
    // bool legal_1 = (_partSize[1] > _min_cell_num && _partSize[1] < _max_cell_num);
    // int max_id = 0;
    // if(legal_0 && legal_1)
    // {
    //     if(_bl[0]->getMaxGainValue() > _bl[1]->getMaxGainValue())
    //     {
    //         max_id = _bl[0]->GetMaxGainNode()->getId();
    //     }
    //     else if(_bl[0]->getMaxGainValue() <= _bl[1]->getMaxGainValue())
    //     {
    //         max_id = _bl[1]->GetMaxGainNode()->getId();
    //     }
    //     else
    //     {
    //         printf("cannot get max gain node\n");
    //         return -1;
    //     }
    // }
    // else if(legal_0 && !legal_1)
    // {
    //     max_id = _bl[0]->GetMaxGainNode()->getId();
    // }
    // else if(!legal_0 && legal_1)
    // {
    //     max_id = _bl[1]->GetMaxGainNode()->getId();
    // }
    // else
    // {
    //     printf("error: max gain node is null\n");
    //     return -1;
    // }
    // return max_id;
    int max_id;
    if(_partSize[1] == _min_cell_num)
        max_id = _bl[0]->GetMaxGainNode()->getId();
    else if(_partSize[0] == _min_cell_num)
        max_id = _bl[1]->GetMaxGainNode()->getId();
    else if(_bl[1]->getMaxGainValue() > _bl[0]->getMaxGainValue())
        max_id = _bl[1]->GetMaxGainNode()->getId();
    else 
        max_id = _bl[0]->GetMaxGainNode()->getId();
    
    
    
    return max_id;
}


void Partitioner::initPartition()
{
    //initialize the partition
    //cell 
    // for (int i = 0; i < _cellArray.size(); ++i) {
    //     // assign the cell to partition A
    //     if(i%2 == 0)
    //     {
    //         _cellArray[i]->setPart(0);
    //         _partSize[0]++;
    //     }
    //     // assign the cell to partition B
    //     else
    //     {
    //         _cellArray[i]->setPart(1);
    //         _partSize[1]++;
    //     }
    // }

    for(int i = 0; i < _cellArray.size()/2; ++i)
    {
        _partSize[0]++;
        _cellArray[i]->setPart(0);
        
    }
    for(int i = _cellArray.size()/2; i < _cellArray.size(); ++i)
    {
        _partSize[1]++;
        _cellArray[i]->setPart(1);
        
    }

    //net
    for(int i = 0; i < _netArray.size(); ++i)
    {
        int part_0_count = 0;
        int part_1_count = 0;
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            if (_cellArray[cellList[j]]->getPart() == 0) {
                part_0_count++;
            }
            else {
                part_1_count++;
            }
        }
        _netArray[i]->setPartCount(0, part_0_count);
        _netArray[i]->setPartCount(1, part_1_count);
    }

}


void Partitioner::initCellGain()
{
    //F(n) = 1 => g(i)++
    //T(n) = 0 => g(i)--
    
    resetBucketlist();
    for(int i = 0; i < _cellArray.size(); ++i)
    {
        Cell* cell = _cellArray[i];
        vector<int> netList = cell->getNetList();
        cell->setGain(0);
        bool from = cell->getPart();
        bool to = !from;
        //iterate through each net
        for(int j = 0; j < netList.size(); ++j)
        {
            Net *net = _netArray[netList[j]];
            if(net->getPartCount(from) == 1) cell->incGain();
            if(net->getPartCount(to) == 0) cell->decGain();
        }

        
        //insert to bucket list corresponding to the gain
        _bl[from]->insertNode(cell->getNode(), cell->getGain());
        
    }
    return;
}

void Partitioner::resetBucketlist()
{
    //reset the bucket list
    for(int i = 0; i < _cellArray.size(); ++i)
    {
        Cell* cell = _cellArray[i];
        cell->setGain(0);
        cell->getNode()->setPrev(nullptr);
        cell->getNode()->setNext(nullptr);
    }

    //bucketlist
    _bl[0]->clearHeadList();
    _bl[1]->clearHeadList();
    return;
}



void Partitioner::cutSize()
{
    //calculate the cut size
    _cutSize = 0;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        if (_netArray[i]->getPartCount(0) > 0 && _netArray[i]->getPartCount(1) > 0) {
            _cutSize++;
        }
    }
    return;
}

void Partitioner::unlockAllCell()
{
    //unlock all the cells
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        _cellArray[i]->unlock();
    }
    return;
}


void Partitioner::updateCellGain(int cellId)
{
    vector<int> netlist = _cellArray[cellId]->getNetList();
    bool from = _cellArray[cellId]->getPart();
    bool to = !from;
    Cell* cell = _cellArray[cellId];
    
    for(int i = 0; i < netlist.size(); ++i)
    {
        Net* net = _netArray[netlist[i]];
        Cell* cell;
        int net_cell_list_size = net->getCellList().size();
        vector<int> cellList = net->getCellList();
        int T = net->getPartCount(to);
        int F = net->getPartCount(from);
        if(T == 0)
        {
            for(int j = 0; j < net_cell_list_size; ++j)
            {
                cell = _cellArray[cellList[j]];
                if(cell->getLock() == false)
                {
                    _cellArray[cellList[j]]->incGain();
                    bucketlistUpdateGain(from, cell);
                }
            }
        }
        else if(T == 1)
        {
            for(int j = 0; j < net_cell_list_size; ++j)
            {
                cell = _cellArray[cellList[j]];
                if(cell->getPart() == to && cell->getLock() == false)
                {
                    cell->decGain();
                    bucketlistUpdateGain(to, cell);
                    break;
                }
            }
        }

        F--;
        T++;


        if(F == 0)
        {
            for(int j = 0; j < net_cell_list_size; ++j)
            {
                cell = _cellArray[cellList[j]];
                if(cell->getLock() == false)
                {
                   
                    cell->decGain();
                    bucketlistUpdateGain(to, cell);
                }
            }
        }
        else if(F == 1)
        {
            for(int j = 0; j < net_cell_list_size; ++j)
            {
                cell = _cellArray[cellList[j]];
                if(cell->getPart() == from && cell->getLock() == false)
                {
                    cell->incGain();
                    bucketlistUpdateGain(from, cell);
                    break;
                }
            }
        }
    }   



}


void Partitioner::bucketlistUpdateGain(bool part, Cell* cell)
{
    _bl[part]->deleteNode(cell->getNode());
    _bl[part]->insertNode(cell->getNode(), cell->getGain());
    return;
}


void Partitioner::moveCell(int cellId)
{
    Cell* cell = _cellArray[cellId];
    bool from = cell->getPart();
    bool to = !from;
    //move the cell
    //partitioner
    _partSize[from]--;
    _partSize[to]++;

    //net
    vector<int> netList = cell->getNetList();
    for(int i = 0; i < netList.size(); ++i)
    {
        Net* net = _netArray[netList[i]];
        net->decPartCount(from);
        net->incPartCount(to);
        
    }
    //cell
    cell->move();

    return;
}





void Partitioner::printSummary() const
{

    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "minimum cell number in each partition: " << _min_cell_num << endl;
    cout << "maximum cell number in each partition: " << _max_cell_num << endl;
    cout << " Balance factor: " << _bFactor << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
