import os
import sys
from command import manager
from utils import client, datainfo
from analyse import analyse

if __name__ == '__main__':
    args = client.parser.parse_args()
  
    # args are required
    if not any(args.__dict__.values()) is True:
        print("Please specify arguments!")
        print("... help message with -h")

    if not os.path.exists('plot') and args.__dict__["which"] == "plot":
        os.makedirs('plot')


    mngr = manager.Run(task = args.__dict__)

    mngr.parse()
    
  
