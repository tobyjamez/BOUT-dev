#!/usr/bin/env python3

  # // N4: matrix(
  # // 	    [0,1,a,a**2/2],
  # // 	    [1,b,b**2/2,b**3/6],
  # // 	    [1,c,c**2/2,c**3/6],
  # // 	    [1,d,d**2/2,d**3/6]
  # // 	    );
  # // iN4:invert(N4);
  # // r:matrix(
  # // 	  [1],
  # // 	  [0],
  # // 	  [0],
  # // 	  [0]
  # // 	  );

def pow(a,b):
    if b == 0:
        return "1"
    if b == 1:
        return a
    else:
        return "%s**%d"%(a,b)
def fac(a):
    if a==0 or a==1:
        return 1
    else:
        assert(a>0)
        return a*fac(a-1)

def dirichlet(order):
    print("A:matrix(")
    for i in range(order):
        print("\t[",end="")
        for j in range(order):
            if j == order-1:
                if i == order-1:
                    end="]\n"
                else:
                    end="],\n"
            else:
                end=", "
            print("%s/%d"%(pow("x%d"%i,j),fac(j)),end=end)
    print(")$")
    print("iA: invert(A)$")
    print("r:matrix([1],")
    for i in range(1,order):
        if i == order-1:
            end="\n"
        else:
            end=",\n"
        print("\t[0]",end=end)
    print(")$")
    print("s:ratsimp(iA.r)$")
    print("""
load (f90) $
""")
    for i in range(order):
        print("f90('fac%d=s[%d,1]);"%(i,i+1))
dirichlet(4)
