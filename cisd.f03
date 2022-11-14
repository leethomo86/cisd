      program truncatedCI 
!
!     This program perfoms a truncated CI calculation.
!
!     L. M. Thompson, 2021
!
      use mqc_gaussian
      use iso_fortran_env
!
!****x* Main/CISD
!*    NAME
!*      Truncated Configuration Interaction Calculator
!*
!*    SYNOPSIS
!*      Computes the energy from a truncated configuration interaction wavefunction. 
!
      implicit none
      character(len=:),allocatable::command,fileName,pLevel,subs_in,help_path
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo 
      integer(kind=int64)::iOut=6,iPrint=1,ival,i,j,flag
      logical::UHF,newNum=.false.
      type(mqc_pscf_wavefunction)::wavefunction
      type(mqc_molecule_data)::moleculeInfo
      type(mqc_twoERIs)::eris,mo_ERIs
      type(mqc_scalar)::Vnn
      type(mqc_determinant)::determinants
      type(mqc_scf_integral)::mo_core_ham
      type(mqc_matrix)::CI_Hamiltonian
      type(mqc_vector)::final_energy,subs
      character(len=10)::val
      integer(kind=int64),dimension(:),allocatable::isubs
!
!*    USAGE
!*      CISD [-f <matrix_file>] [-p <print_level>] [-s <substitutions>] [--help]
!*
!*    OPTIONS
!* 
!
!     Print program information.
!
      Write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',74),NEW_LINE('a'), &
          repeat(' ',11),'Truncated Configuration Interaction Energy Calculator',NEW_LINE('a'), &
          ' ',repeat('*',74),NEW_LINE('a'), &
          NEW_LINE('a'),repeat(' ',30),'Version 22.11.14',NEW_LINE('a'),NEW_LINE('a'),&
          ' L. M. Thompson, Louisville KY, 2021.',NEW_LINE('a')
!
!     Parse input options.
!
      j = 1
      do i=1,command_argument_count()
        if(i.ne.j) cycle
        call mqc_get_command_argument(i,command)
        if(command.eq.'-f') then
!
!*      -f matrix_file                   Input matrix file with initial set of molecular orbitals.
!*
          call mqc_get_command_argument(i+1,fileName)
          j = i+2
        elseif(command.eq.'-p') then
!
!*      -p print_level                   Verbosity of output. Default print level is 1. Options
!*                                       0-4.
!*
          call mqc_get_command_argument(i+1,pLevel)
          read(pLevel,'(I1)') iPrint
          j = i + 2
        elseIf(command.eq.'-s') then
!
!*      -s substitutions                 Substitution levels permitted in truncated CI calculation.
!*                                       Input is a comma delimited list of numbers inside square 
!*                                       brackets. The reference is always included and any level
!*                                       and combination of truncations can be specified. 
!*
!*                                       E.g. [1,2,4] specifies truncated CI with all single, double 
!*                                       and quadruple substitutions.
!*
          call mqc_get_command_argument(i+1,subs_in)
          j = i+2
        elseIf(command.eq.'--help') then
!
!*      --help                           Output help documentation to terminal.
!*
          if(command_argument_count().gt.1) call mqc_error_I('Help output requested with multiple arguments',6, &
            'command_argument_count()',command_argument_count())
          call mqc_get_command_argument(0,help_path)
          help_path = 'less ' // trim(help_path(1:scan(help_path,'/',.true.))) // 'doc/CISD.txt'
          call execute_command_line(help_path,exitstat=flag)
          if(flag.ne.0) call mqc_error('Help output command failed')
          stop
        else
          call mqc_error_A('Unrecognised input flag',6,'command',command)
        endIf
      endDo
!
      if(.not.allocated(subs_in)) then
        write(iOut,'(A)') 'Defaulting to CISD.'
        subs = [1,2]
      else
        do i = 1,len(subs_in)
          select case (subs_in(i:i))
          case('[','\(')
            if(i.eq.1) then
              newNum = .true.
              cycle
            else
              call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
            endIf
          case(']','\)')
            if(i.eq.len(subs_in).and..not.newNum.and.i.ne.1) then
              read(val,'(I10)') ival
              call subs%push(ival)
              newNum = .true.
              cycle
            else
              call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
            endIf
          case(',',' ')
            if(i.eq.1.or.i.eq.len(subs_in).or.i.eq.2.or.i.eq.len(subs_in)-1.or.newNum) then
              call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
            else
              read(val,'(I10)') ival
              call subs%push(ival)
              newNum = .true.
              cycle
            endIf
          case('0':'9')
            if(i.eq.1.or.i.eq.len(subs_in)) then
              call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
            else
              if(newNum) then
                val = subs_in(i:i)
              else
                val = trim(val)//subs_in(i:i)
              endIf
              newNum = .false.
              cycle
            endIf
          case default
            call mqc_error_A('Substitution level input format incorrect',6,'subs_in',subs_in)
          end select
        endDo
      endIf
!
      call subs%print(6,'Permitted substitution levels')
!
      call fileInfo%load(filename)
      call fileInfo%getMolData(moleculeInfo)
      call moleculeInfo%print(iOut)
      call fileInfo%getESTObj('wavefunction',wavefunction)
      if(iPrint.ge.2) call wavefunction%print(iOut,'all')
      call fileInfo%get2ERIs('regular',eris)
      if(iPrint.ge.4) call eris%print(iOut,'AO 2ERIs')
!
      if(wavefunction%wf_type.eq.'U') then
        UHF = .true.
        write(iOut,*) 'found UHF wavefunction'
      elseIf(wavefunction%wf_type.eq.'R') then
        UHF = .False.
        write(iOut,*) 'found RHF wavefunction'
      else
        call mqc_error_A('Unsupported wavefunction type in fullci',iOut, &
          'Wavefunction%wf_type',wavefunction%wf_type)
      endIf 
!
      if (wavefunction%wf_complex) call mqc_error('Complex wavefunctions unsupported in fullci')
!
!     Compute the nuclear-nuclear repulsion energy.
!
      Vnn = mqc_get_nuclear_repulsion(moleculeInfo)
      call Vnn%print(iOut,'Nuclear Repulsion Energy (au)') 
!
!     Generate Slater Determinants       
!
      if(iPrint.ge.2) write(iOut,*) 'Building Determinant Strings'
      isubs = [(i, i=1,int(maxval(subs)))]
      call trci_dets_string(iOut,iPrint,wavefunction%nBasis,wavefunction%nAlpha, &
        Wavefunction%nBeta,isubs,determinants)
!
!     Transform one and two-electron integrals to MO basis
!
      if(iPrint.ge.2) write(iOut,*) 'Transforming MO integrals'
      mo_core_ham = matmul(transpose(wavefunction%MO_Coefficients),matmul(wavefunction%core_Hamiltonian, &
          Wavefunction%MO_Coefficients))
      if(IPrint.ge.3) call mo_core_ham%print(iOut,'MO Basis Core Hamiltonian') 
      call twoERI_trans(iOut,iPrint,wavefunction%MO_Coefficients,ERIs,mo_ERIs)
!
!     Generate Hamiltonian Matrix
!
      if(iPrint.eq.2) write(iOut,*) 'Building CI Hamiltonian'
      call subs%unshift(0)
      isubs = subs
      call mqc_build_ci_hamiltonian(iOut,iPrint,wavefunction%nBasis,determinants, &
        mo_core_ham,mo_ERIs,UHF,CI_Hamiltonian,isubs)
!
!     Diagonalize Hamiltonian
!
      if(iPrint.ge.2) write(iOut,*) 'Diagonalizing CI Hamiltonian'
      call CI_Hamiltonian%diag(wavefunction%pscf_energies,wavefunction%pscf_amplitudes)
      if(iPrint.ge.3) then 
        call wavefunction%pscf_amplitudes%print(iOut,'CI Eigenvectors')
        call wavefunction%pscf_energies%print(iOut,'CI Eigenvalues')
      endIf
!
      final_energy = Vnn + wavefunction%pscf_energies
      if(iPrint.ge.2) call final_energy%print(iOut,'Final Energy') 
      call mqc_print(final_energy%at(1),iOut,'Ground State Energy',.true.,.true.,'F16.10')
!
!*    NOTES
!*      Compilation of this program requires the MQC library (https://github.com/MQCPack/mqcPack)
!*      and the gauopen utility (http://gaussian.com/g16/gauopen.zip) and compilation with the
!*      f08 standard.
!*
!*      Compilation tested using: gfortran 12.2.0
!*
!*      Note that subroutine Wr_LCBuf needs modifying in gauopen/qcmatrix.F as follows:
!*        line 58:       LenBX = (LenBuf/(2*abs(NR)))*abs(NR)
!*        line 60:       Call Wr_CBuf(IU,NTot*abs(NR),LenBX,X)
!*
!*      Documentation generated with robodoc. To update documentation edit robodoc.rc to
!*      determine documentation output type and then run robodoc at the command line in the
!*      main directory.
!*
!*    AUTHORS
!*      Lee M. Thompson, University of Louisville, lee.thompson.1@lousiville.edu
!*
!*    COPYRIGHT
!*      (c) 2021 by Lee M. Thompson distributed under terms of the MIT license.
!*
!****
!
 999  End Program truncatedCI
