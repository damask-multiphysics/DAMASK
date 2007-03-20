
!##############################################################
      MODULE mesh     
!##############################################################

      use prec, only: pRe,pIn
      implicit none

!     ---------------------------
!	_Nelems	: total number of elements in mesh
!	_Nnodes	: total number of nodes in mesh
!	_maxNnodes	: max number of nodes in any element
!	_maxNips	: max number of IPs in any element
!	_element	: type, material, node indices
!	_node		: x,y,z coordinates (initially!)
!	_nodeIndex	: count of elements containing node,
!			  [element_num, node_index], ...
!	_envIP	: 6 neighboring IPs as [element_num, IP_index]
!			  order is +x, +y,+z, -x, -y, -z in local coord
!     ---------------------------
      integer(pIn) mesh_Nelems, mesh_Nnodes, mesh_maxNnodes,mesh_maxNips
      integer(pIn), allocatable :: mesh_element (:,:)
      integer(pIn), allocatable :: mesh_nodeIndex (:,:)
      integer(pIn), allocatable :: mesh_envIP (:,:)
      real(pRe), allocatable :: mesh_node (:,:)

	CONTAINS
!     ---------------------------
!	subroutine mesh_init()
!	subroutine mesh_parse_inputFile()
!     ---------------------------

	
!     ***********************************************************
!     initialization 
!     ***********************************************************
	SUBROUTINE mesh_init ()

	mesh_Nelems = 0_pIn
	mesh_Nnodes = 0_pIn
	mesh_maxNips = 0_pIn
	mesh_maxNnodes = 0_pIn
	call mesh_parse_inputFile ()

	END SUBROUTINE
	
!     ***********************************************************
!     parsing of input file 
!     ***********************************************************
	FUNCTION mesh_parse_inputFile ()

      use prec, only: pRe,pIn
      use IO
      
      implicit none

	logical mesh_parse_inputFile
      integer(pIn) i,j,positions(10*2+1)
      integer(pIn) elem_num,elem_type,Nnodes,node_num,num_ip,mat,tp(70,2)

!     Set a format to read the entire line (max. len is 80 characters)
  610 FORMAT(A80)

	if (.not. IO_open_inputFile(600)) then
	  mesh_parse_inputFile = .false.
	  return
	endif

      do while(.true.)
        read(600,610,end=620) line

	  positions = IO_stringPos(line,3)
	  select case (IO_stringValue(line,positions,1)
!-----------------------------------
	  case ('sizing')
!-----------------------------------
	    mesh_Nelems = IO_intValue(line,positions,2)
	    mesh_Nnodes = IO_intValue(line,positions,3)
!-----------------------------------
	  case ('elements')
!-----------------------------------
	    select case (IO_intValue(line,positions,2)) ! elem type
	      case (3) ! 2D Triangle
	        mesh_maxNips = max(3,mesh_maxNips)
	        mesh_maxNnodes = max(3,mesh_maxNnodes)
	      case (6) ! 2D Quad.
	        mesh_maxNips = max(4,mesh_maxNips)
	        mesh_maxNnodes = max(4,mesh_maxNnodes)
	      case (7) ! 3D hexahedral
	        mesh_maxNips = max(8,mesh_maxNips)
	        mesh_maxNnodes = max(8,mesh_maxNnodes)
	      case default
	        mesh_maxNips = max(8,mesh_maxNips)
	        mesh_maxNnodes = max(8,mesh_maxNnodes)
	    end select
!-----------------------------------
	  case ('connectivity')
!-----------------------------------
          allocate (mesh_element(mesh_Nelems,2+mesh_maxNips))
          allocate (mesh_nodeIndex (mesh_Nnodes,1+mesh_maxNnodes*2)
          allocate (mesh_envIP(mesh_Nelems,mesh_maxNips,6,2))
	    mesh_element = 0_pIn
	    mesh_nodeIndex = 0_pIn
	    mesh_envIP = 0_pIn

!	MISSING: setting up of envIP

          read(600,610,end=620) line ! skip line ??
          do i=1,mesh_Nelems
            read(600,610,end=620) line
	      positions = IO_stringPos(line,0) ! find all chunks
	      elem_num  = IO_intValue(line,positions,1)
	      elem_type = IO_intValue(line,positions,2)
	      select case (elem_type)
	      case (3) ! 2D Triangle
	        Nnodes = 3
	      case (6) ! 2D Quad.
	        Nnodes = 4
	      case (7) ! 3D hexahedral
	        Nnodes = 8
	      case default
	        Nnodes = 8
	      end select
            mesh_element(elem_num,1) = elem_type
	      do j=1,Nnodes  ! store all node indices
	        node_num = IO_intValue(line,positions,2+j)
              mesh_element(elem_num,1+j) = node_num
              mesh_nodeIndex(node_num,1) = mesh_nodeIndex(node_num,1)+1 ! inc count
              mesh_nodeIndex(node_num,mesh_nodeIndex(node_num,1)*2  ) = elem_num
              mesh_nodeIndex(node_num,mesh_nodeIndex(node_num,1)*2+1) = j
            end do
          end do
!-----------------------------------
	  case ('coordinates')
!-----------------------------------
          allocate (mesh_node(mesh_Nnodes,3)) ! x,y,z, per node

          read(600,610,end=620) line ! skip line ??
          do i=1,mesh_Nnodes
            read(600,610,end=620) line
	      positions = IO_stringPos(line,0) ! find all (4) chunks
	      node_num  = IO_intValue(line,positions,1)
	      do j=1,3  ! store x,y,z coordinates
              mesh_node(node_num,j) = IO_floatValue(line,positions,1+j)
            end do
          end do
!-----------------------------------
	  case ('hypoelastic')
!-----------------------------------
 
!       ***************************************************
!       Search for key word "hypoelastic".  
!         This section contains the # of materials and
!           the element range of each material
!       ***************************************************
        ELSE IF( line(1:11).eq.'hypoelastic' )THEN
          mat=0
          flag=0
          DO WHILE( line(1:8).ne.'geometry' )
            READ(600,610,END=620) line
            i=1
            DO WHILE( i.le.len(line)-8 )
              IF( line(i:i+2).eq.'mat' )THEN
                mat=mat+1
                flag=1
              END IF
              i=i+1
            END DO
            IF( flag.eq.1 )THEN
              flag=0
              READ(600,610,END=620) line
              READ(600,610,END=620) line
              i=1
              DO WHILE( line(i:i).eq.' ')
                i=i+1
              END DO
              left=i
              DO WHILE( line(i:i).ne.' ')
                i=i+1
              END DO
              right=i-1
              READ(UNIT=line(left:right), FMT=' (I5) ') tp(mat,1)
              DO WHILE( (line(i:i).eq.' ').or.
     &                  (line(i:i).eq.'t').or.
     &                  (line(i:i).eq.'o') )
                    i=i+1
              END DO
              left=i
              DO WHILE( line(i:i).ne.' ')
                i=i+1
              END DO
              right=i-1
              READ(UNIT=line(left:right), FMT=' (I5) ') tp(mat,2)
            END IF
          END DO
          WRITE(6,*) 'mat: ',mat,' ',tp(1,1),' ',tp(1,2)

	  end select
      END DO

!     Code jumps to 620 when it reaches the end of the file
  620 continue

      WRITE(6,*) 'Finished with .dat file'
      CALL FLUSH(6)

	END FUNCTION
	
	END MODULE mesh
	