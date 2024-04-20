import time
import shutil
import os
import subprocess

# ids = [6,7,8]
# for id in ids:
#     #  cmd = f"matlab -batch run('id={id};easycall.m;quit') -logfile DMTlog/testout{id}.log"
#      cmd = f"matlab -nodesktop -nosplash -r \"id={id};easycall;quit\" -logfile DMTlog/DMT_spiral_detection_sub{id:02}.log"
#      print(cmd)
#      os.system(cmd)



# ids = [6, 7, 8]
ids = range(9, 23)
for id in ids:
    cmd = f"matlab -batch \"id={id};easycall;quit\" -logfile DMTlog/DMT_spiral_detection_sub{id:02}.log"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)
