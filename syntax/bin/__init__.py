from syntax import *

        
if __name__=='__main__':
    netlist = FastHenry()
    netlist.append( Title('my script') )
    netlist.append( Units('um') )
    netlist.append( Default(z=0, sigma=5e7) )
    netlist.append( Freq(1e8, 10e9, 1) )
    n1 = Node(0, 0)
    n2 = Node(500, 0)
    n3 = Node(500, 500)
    n4 = Node(0, 500)
    n5 = Node(0, 1)
    n6 = Node(100, 100)
    n7 = Node(400, 100)
    n8 = Node(400, 400)
    n9 = Node(100, 400)
    n10 = Node(100, 101)

    netlist.extend([n1,n2,n3,n4,n5,n6,n7,n8,n9,n10])
    netlist.append( Equiv(n1, n6) )
    netlist.append( Equiv(n5, n10) )

    nhinc = 5
    nwinc = nhinc*3

    width = 10
    thickness = 3.4
    netlist.append( Segment(n1,n2,w=width,h=thickness,nhinc=nhinc,nwinc=nwinc) )
    netlist.append( Segment(n2,n3,w=width,h=thickness,nhinc=nhinc,nwinc=nwinc) )
    netlist.append( Segment(n3,n4,w=width,h=thickness,nhinc=nhinc,nwinc=nwinc) )
    netlist.append( Segment(n4,n5,w=width,h=thickness,nhinc=nhinc,nwinc=nwinc) )
    
    width = 8
    thickness = 3.4
    netlist.append( Segment(n6,n7,w=width,h=thickness,nhinc=nhinc,nwinc=nwinc) )
    netlist.append( Segment(n7,n8,w=width,h=thickness,nhinc=nhinc,nwinc=nwinc) )
    netlist.append( Segment(n8,n9,w=width,h=thickness,nhinc=nhinc,nwinc=nwinc) )
    netlist.append( Segment(n9,n10,w=width,h=thickness,nhinc=nhinc,nwinc=nwinc) )

    netlist.append( Port(n1, n5) )
    #netlist.append( Port(n6, n10) )
    
    netlist.run(verbose=False, refinement=1)
    print netlist.stdout
