mol delete top
mol load vtf vmd-trajectory.vtf
mol delrep 0 top
display resetview
pbc box_draw -color purple -width 3
set typ51 [atomselect top "name 51 "]
$typ51 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 51
mol material Opaque
mol color ColorID 0
mol addrep top
 
set typ20 [atomselect top "name 20 "]
$typ20 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 20
mol material Opaque
mol color ColorID 1
mol addrep top
 
set typ38 [atomselect top "name 38 "]
$typ38 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 38
mol material Opaque
mol color ColorID 2
mol addrep top
 
set typ25 [atomselect top "name 25 "]
$typ25 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 25
mol material Opaque
mol color ColorID 3
mol addrep top
 
set typ40 [atomselect top "name 40 "]
$typ40 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 40
mol material Opaque
mol color ColorID 4
mol addrep top
 
set typ37 [atomselect top "name 37 "]
$typ37 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 37
mol material Opaque
mol color ColorID 5
mol addrep top
 
set typ45 [atomselect top "name 45 "]
$typ45 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 45
mol material Opaque
mol color ColorID 6
mol addrep top
 
set typ41 [atomselect top "name 41 "]
$typ41 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 41
mol material Opaque
mol color ColorID 7
mol addrep top
 
set typ49 [atomselect top "name 49 "]
$typ49 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 49
mol material Opaque
mol color ColorID 8
mol addrep top
 
set typ18 [atomselect top "name 18 "]
$typ18 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 18
mol material Opaque
mol color ColorID 9
mol addrep top
 
set typ19 [atomselect top "name 19 "]
$typ19 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 19
mol material Opaque
mol color ColorID 10
mol addrep top
 
set typ50 [atomselect top "name 50 "]
$typ50 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 50
mol material Opaque
mol color ColorID 11
mol addrep top
 
set typ24 [atomselect top "name 24 "]
$typ24 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 24
mol material Opaque
mol color ColorID 12
mol addrep top
 
set typ48 [atomselect top "name 48 "]
$typ48 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 48
mol material Opaque
mol color ColorID 13
mol addrep top
 
set typ36 [atomselect top "name 36 "]
$typ36 set radius 1
mol representation VDW 1.000000 16.000000
mol selection name 36
mol material Opaque
mol color ColorID 14
mol addrep top
 
animate goto 0
color Display Background white
axes location off
