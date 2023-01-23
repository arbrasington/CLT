# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:00:15 2023

@author: ALEXRB
"""


class Node:
    def __init__(self, key_dict):
        del key_dict['__class__']
        del key_dict['self']

        for key, value in key_dict.items():
            try:
                setattr(self, key, value)
            except:
                print('Ivalid keyword argument')
