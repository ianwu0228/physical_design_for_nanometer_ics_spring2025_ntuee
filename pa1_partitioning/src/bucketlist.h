//bucketlist.h

#ifndef BUCKETLIST_H
#define BUCKETLIST_H
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "cell.h"
using namespace std;

class Bucketlist
{

    public:
        //constructor 
        Bucketlist(int max_pin_num)
        {   
            max_pin_num_ = max_pin_num;
            max_gain_ = -max_pin_num_;
            for(int i = -max_pin_num_; i <= max_pin_num_; ++i)
            {
                head_list_.insert(pair<int, Node*>(i, new Node()));
            }
        }

        //destructor
        ~Bucketlist()
        {
            for(int i = -max_pin_num_; i <= max_pin_num_; ++i)
            {
                delete head_list_[i];
            }
        }


        void insertNode(Node* node, int gain)
        {
            
            
            if(gain > max_gain_)
            {
                max_gain_ = gain;
            }

            if(head_list_[gain]->getNext() == nullptr)
            {
                head_list_[gain]->setNext(node);
                node->setPrev(head_list_[gain]);
                node->setNext(nullptr);
            }
            else
            {
                Node* temp = head_list_[gain]->getNext();
                head_list_[gain]->setNext(node);
                node->setPrev(head_list_[gain]);
                node->setNext(temp);
                temp->setPrev(node);
            }
            return;
            
        }

        void deleteNode(Node* node)
        {
            
            if(node == nullptr)
            {
                printf("node deleted is null, error\n");
                return;

            }
            //check if the deleted node is the only max gain node left
            if(head_list_[max_gain_]->getNext() == node && node->getNext() == nullptr)
            {
                // printf("deleting cell id = %d\n", node->getId());   
                // printf("delete head node\n");
                Node* max_gain_node = head_list_[max_gain_]->getNext();
                head_list_[max_gain_]->setNext(nullptr);
                max_gain_node->setPrev(nullptr);
                max_gain_node->setNext(nullptr);
                while(max_gain_ > -max_pin_num_ && head_list_[max_gain_]->getNext() == nullptr)
                {
                    max_gain_--;
                    // printf("-- max gain updating = %d\n", max_gain_);
                }
                return;
            }

            //else we can treat all the cases the same way since the max gain will not be altered
            //node is the last element in the list
            if(node->getNext() == nullptr)
            {
                // if(node->getPrev() == nullptr)printf("error: node is null\n");
                Node* prev_node = node->getPrev();
                if(prev_node == nullptr)
                {
                    printf("error: node is null0101\n");
                    printf("cell id = %d\n", node->getId());
                }
               
                prev_node->setNext(nullptr);
                node->setPrev(nullptr);
                node->setNext(nullptr);
            }
            //node is not the last element in the list
            else
            {
                Node* temp = node->getNext();
                node->getPrev()->setNext(temp);
                temp->setPrev(node->getPrev());
                node->setPrev(nullptr);
                node->setNext(nullptr);
            }
            
        }


        // void deleteMaxGainNode()
        // {
        //     if(head_list_[max_gain_]->getNext()->getNext() == nullptr)//only max gain node
        //     {
        //         Node* max_gain_node = head_list_[max_gain_]->getNext();
        //         head_list_[max_gain_]->setNext(nullptr);
        //         max_gain_node->setPrev(nullptr);
        //         max_gain_node->setNext(nullptr);
        //         while(max_gain_ >= -max_pin_num_ && head_list_[max_gain_]->getNext() == nullptr)
        //         {
        //             max_gain_--;
        //         }
        //     }
        //     //not the only one for the same max gain
        //     else
        //     {
        //         Node* max_gain_node = head_list_[max_gain_]->getNext();
        //         head_list_[max_gain_]->setNext(max_gain_node->getNext());
        //         max_gain_node->getNext()->setPrev(head_list_[max_gain_]);
        //         max_gain_node->setPrev(nullptr);
        //         max_gain_node->setNext(nullptr);
        //     }
        // }


        Node* GetMaxGainNode()
        {
            Node* node;
            if(head_list_[max_gain_]->getNext() == nullptr)
            {
                printf("erroe: max gain node is null\n");
                node = nullptr;
                
            }
            else
            {
                node = head_list_[max_gain_]->getNext();
            }
            return node;
        }

        void printBucketlist()
        {
            printf("max gain = %d\n", max_gain_);
            for(int i = -max_pin_num_; i <= max_pin_num_; ++i)
            {
                printf("gain = %d: ", i);
                Node* node = head_list_[i]->getNext();
                while(node != nullptr)
                {
                    printf("%d -> ", node->getId());
                    node = node->getNext();
                }
                printf("\n");
            }
        }

        int getMaxGainValue()
        {
            return max_gain_;
        }


        void clearHeadList()
        {
            for(int i = -max_pin_num_; i <= max_pin_num_; ++i)
            {
                head_list_[i]->setNext(nullptr);
                head_list_[i]->setPrev(nullptr);
            }
        }

        

        


    private:
        int max_gain_;
        int max_pin_num_;
        map<int, Node*> head_list_;




};



#endif
