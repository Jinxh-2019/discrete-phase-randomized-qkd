import sys
import time
def progress_bar(func):
    for i in range(1,101):
        print("\r",end='')
        print("Downloading progress: {}%",'▋'*(i//2) ,end='')
        sys.stdout.flush()
        func()
