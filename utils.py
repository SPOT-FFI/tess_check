"""
Module to contain miscellaneous utility functions
"""


def list_to_string(list,delim=' '):
  list_as_string = list[0]
  for i in range(len(list)-1):
    id = i+1
    list_as_string += (delim + list[id])
  return list_as_string
