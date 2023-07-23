program prepare
implicit none

integer                 :: i, j, k, l, max
integer, parameter      :: max_seq = seq_number         ! Total number of sequences
integer, parameter      :: max_hb = hb_number           ! Total number of HBonds
integer, parameter      :: l_seq = 3                    ! Sequence ranges
character*3             :: res_seq(max_seq, 2)          ! Residues for sequences (1:origin/2:amber)
integer                 :: n_seq(max_seq, 2)            ! Numbering of sequences (1:origin/2:amber)
character*1             :: c_seq(max_seq, 2)            ! Chain number of sequences (1:origin/2:amber)
character*1             :: ss_seq(max_seq)              ! Secondary structure of sequences
character*3             :: res(2)                       ! Residues for HBonds (1:donor/2:acceptor)
integer                 :: n(2)                         ! (I)Numbering of HBonds residues (1:donor/2:acceptor)
integer                 :: m(2)                         ! (O)Numbering of HBonds residues (1:donor/2:acceptor)
character*4             :: at(2)                        ! Atom type of HBonds residues (1:donor/2:acceptor)
character*1             :: a(2)                         ! Atom of HBonds residues (1:donor/2:acceptor)
double precision        :: r                            ! HBond length
integer                 :: c(2)                         ! Charges of HBonds (1:donor/2:acceptor)
character*1             :: p(2)                         ! Postions of HBonds (1:donor/2:acceptor)
character*1             :: ss(2)                        ! Secondary structure of HBonds (1:donor/2:acceptor)
character*3             :: seq(-l_seq : l_seq, 2)       ! Sequences (-3/+3) of HBonds (1:donor/2:acceptor)
integer                 :: error
character*99            :: fmt = "(A3,X,A1,X,A3,X,A1,2I3,X,A1,X,A1,X,A1,X,A1,14A4,F7.2,X,A4,X,A4,2I5)"

open (11, file = 'seq_origin')
open (12, file = 'seq_amber')
open (13, file = 'temp_ss')
open (21, file = 'hbond')
open (99, file = 'output.log')

error = 0
!!!!! Read the sequences !!!!!
do i = 1, max_seq
   read (11, "(A3, X, A1, I4)") res_seq(i, 1), c_seq(i, 1), n_seq(i, 1)
   read (12, "(A3, X, A1, I4)") res_seq(i, 2), c_seq(i, 2), n_seq(i, 2)
   if (res_seq(i, 1) /= res_seq(i, 1)) then
      error = 1
   end if
   if (i /= 1) then
      if (c_seq(i, 1) == c_seq(i - 1, 1) .and. n_seq(i, 1) == n_seq(i - 1, 1)) then
         error = 2
      end if
   end if
end do
close (11)
close (12)

if (error /= 0) then
   open (90, file = 'error')
   if (error == 1) then
      write (90, *) 'Error 1: Non-Amber-compatible residues marked as ATOM. Please change it into HETATM.'
   end if
   if (error == 2) then
      write (90, *) 'Error 2: Please make sure one residue sequence number only contains one type of residue.'
   end if
   close (90)
else
!!!!! Read the secondary structures !!!!!
max = int(max_seq / 50)
do i = 1, max
   read (13, *)
   read (13, "(50A1)") (ss_seq((i - 1) * 50 + k), k = 1, 50)
   read (13, *)
end do
if (max_seq - max * 50 /= 0) then
read (13, *)
read (13, "(50A1)") (ss_seq(max * 50 + k), k = 1, max_seq - max * 50)
read (13, *)
end if
do i = 1, max_seq
   if (ss_seq(i) == ' ') then
      ss_seq(i) = 'C'
   end if
end do
close (13)

!!!!! Read the HBonds and make make the lines of parameters !!!!!
do i = 1, max_hb
   read (21, *) res(2), n(2), at(2), res(1), n(1), at(1), r
   do k = 1, 2
      a(k) = at(k)(1:1)
      if (at(k)(2:3) == '  ' .or. at(k)(2:3) == 'XT') then
         p(k) = 'B'
         c(k) = 0
      else
         p(k) = 'S'
         c(k) = 0
         if (res(k) == 'ARG' .or. res(k) == 'HIP' .or. res(k) == 'LYS') then
            c(k) = 1
         end if
         if (res(k) == 'ASP' .or. res(k) == 'GLU' .or. res(k) == 'CYM' .or. res(k) == 'TYM') then
            c(k) = -1
         end if
      end if
      if (res(k) == 'HIP' .or. res(k) == 'HIE' .or. res(k) == 'HID') then
         res(k) = 'HIS'
      end if
      if (res(k) == 'CYX') then
         res(k) = 'CYS'
      end if
      seq(0, k) = res(k)
      do j = 1, max_seq
         if (n(k) == n_seq(j, 2)) then
            ss(k) = ss_seq(j)
            m(k) = n_seq(j, 1)
            do l = 1, l_seq
               if (j - l > 0) then
                  if (n_seq(j, 1) - n_seq(j - l, 1) == l .and. c_seq(j, 1) == c_seq(j - l, 1)) then
                     seq(0 - l, k) = res_seq(j - l, 2)
                  else
                     seq(0 - l, k) = 'XXX'
                  end if
               else
                  seq(0 - l, k) = 'XXX'
               end if
               if (j + l <= max_seq) then
                  if (n_seq(j + l, 1) - n_seq(j, 1) == l .and. c_seq(j + l, 1) == c_seq(j, 1)) then
                     seq(0 + l, k) = res_seq(j + l, 2)
                  else
                     seq(0 + l, k) = 'XXX'
                  end if
               else
                  seq(0 + l, k) = 'XXX'
               end if
            end do
         end if
      end do
   end do
   k = 0
   do l = -l_seq, l_seq
      if (seq(l, 1) == 'XXX' .or. seq(l, 2) == 'XXX') then
         k = 1
      end if
   end do
   if (k == 0) then
if (p(1) == 'S') then
write (99, fmt) res(1),a(1),res(2),a(2),c(1:2),p(1:2),ss(1:2),seq(-l_seq:l_seq,1),seq(-l_seq:l_seq,2),r,at(1:2),m(1:2)
end if
   end if
end do
end if
close (21)
close (99)

end program prepare
