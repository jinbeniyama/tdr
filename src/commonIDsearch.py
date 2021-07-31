#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Search common IDs.
"""
from argparse import ArgumentParser as ap

def main(args=None):
  """This is the main function called by the `commonIDsearch` script.

  Parameters
  ----------
  args : argparse.Namespace
    Arguments passed from the command-line as defined below.
  """
  parser = ap(description="Convert ID to fitslist")
  parser.add_argument(
    "ID", nargs="*", 
    help="ID text")
  parser.add_argument(
    "--txt", nargs="*", 
    help="original text")
  parser.add_argument(
    "--pre", type=str, default=None, 
    help="word before ID")
  parser.add_argument(
    "--post", type=str, default=None, 
    help="word after ID")
  args = parser.parse_args()
  
  # Save IDs
  ID_list_all = []
  for ID in args.ID:
    ID_list = []
    with open(ID, "r") as f:
      f = f.readlines()
      for line in f:
        ID = line.strip("\n")
        ID_list.append(ID)
      ID_list_all.append(ID_list)
  
  # Extract common ID
  for n in range(len(ID_list_all)-1):
    ID_common = list(set(ID_list_all[n]) & (set(ID_list_all[n+1])))
  
  # Sort common ID
  ID_common.sort()
  print(f"common ID N={len(ID_common)}")

  # Create common fits list using common ID
  # output {pre} + {commonID} + {post} 
  for ID in ID_common:
    print(f"{args.pre}{ID}{args.post}")


if __name__ == "__main__":
  parser = ap(description="Convert ID to fitslist")
  parser.add_argument(
    "ID", nargs="*", 
    help="ID text")
  parser.add_argument(
    "--txt", nargs="*", 
    help="original text")
  parser.add_argument(
    "--pre", type=str, default=None, 
    help="word before ID")
  parser.add_argument(
    "--post", type=str, default=None, 
    help="word after ID")
  args = parser.parse_args()
  
  main(args)
