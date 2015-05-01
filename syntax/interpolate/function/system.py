import re, time, os
from subprocess import Popen, PIPE, STDOUT

def find(path='.', regex='*', ctime=0):
    r = []
    regex = str(regex).strip()
    if regex == '*': regex = ''
    now = time.time()
    for filename in os.listdir(path):
        try:
            if re.search(regex, filename):
                tmtime = os.path.getmtime(os.path.join(path, filename))
                if ctime>0 and int((now-tmtime)/3600/24) > ctime:
                    r.append(os.path.join(path, filename))
                elif ctime<0 and int((now-tmtime)/3600/24) < ctime:
                    r.append(os.path.join(path, filename))
                elif ctime==0:
                    r.append(os.path.join(path, filename))
        except:
            pass
    return r


def removedirs(*files):
    for i, file in enumerate(files):
        try:
            os.system('/bin/rm -rf %s > /dev/null 2>&1'%file)
	except:
	    pass
#
#        try:
#            os.removedirs(file)
#            continue
#        except:
#            pass
#        try:
#            os.remove(file)
#            continue
#        except:
#            pass


def source(filename):
        cmd = "source {filename}; env".format(filename=filename)
        p = Popen(cmd, executable='/bin/tcsh', stdout=PIPE, stderr=STDOUT, shell=True, env=os.environ)
        stdout = p.communicate()[0].splitlines()
        for line in stdout:
            if re.search('[0-9a-zA-Z_-]+=\S+', line):
                key, value = line.split("=", 1)
                os.environ[key] = value
