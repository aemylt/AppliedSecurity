import sys, subprocess

def interact( G ) :
  # Send      G      to   attack target.
  target_in.write( "%s\n" % ( G ) ) ; target_in.flush()

  # Receive ( t, r ) from attack target.
  t = int( target_out.readline().strip() )
  r = int( target_out.readline().strip() )

  return ( t, r )

def attack() :
  # Select a hard-coded guess ...
  G = list("a")
  r = 0
  old_t = 0
  i = 0

  while r != 1:
    i += 1
    # ... then interact with the attack target.
    ( t, r ) = interact( "".join(G) )
    if t == 0:
      G.append("a")
    elif t != old_t:
      old_t = t
    else:
      G[t - 1] = chr(ord(G[t - 1]) + 1)

  # Print all of the inputs and outputs.
  print "G = %s" % ( "".join(G) )
  print "t = %d" % ( t )
  print "r = %d" % ( r )
  print "i = %d" % ( i )

if ( __name__ == "__main__" ) :
  # Produce a sub-process representing the attack target.
  target = subprocess.Popen( args   = sys.argv[ 1 ],
                             stdout = subprocess.PIPE, 
                             stdin  = subprocess.PIPE )

  # Construct handles to attack target standard input and output.
  target_out = target.stdout
  target_in  = target.stdin

  # Execute a function representing the attacker.
  attack()
