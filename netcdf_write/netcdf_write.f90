
module netcdf_write
  use netcdf
  
  implicit none
  private

  integer :: fid=-1, xid=-1, yid=-1, rdimid=-1, timeid
  type vartype
     character (len=256) :: nm
     integer :: ID
  end type vartype

  type(vartype) :: var(50)
  integer :: nvars=0

  public :: write_nc

  interface write_nc
     module procedure write_nc_2d
     module procedure write_nc_2d_1 
  end interface write_nc
  
contains

  subroutine write_nc_2d(xaxis, yaxis, field, name, ntime)
     real, intent(in) :: xaxis(:), yaxis(:), field(:,:)
     character (len=*), intent(in) :: name
     integer, optional :: ntime
      
     !local vars
     integer :: xvarid, yvarid, stat, varid
     integer :: ntime1, i
      
     ntime1 = 1
     if (present(ntime)) ntime1=ntime
      
     if (fid<0) then
        stat = nf90_create('output.nc', nf90_clobber, fid)
        call handle_err(stat)
        stat = nf90_def_dim(fid, "X", size(xaxis), xid)
        call handle_err(stat)
        stat = nf90_def_dim(fid, "Y", size(yaxis), yid)
        call handle_err(stat)
        stat = nf90_def_dim(fid, "time", nf90_unlimited, rdimid)
        call handle_err(stat)
         
        stat = nf90_def_var(fid, "X", nf90_float, xid, xvarid)
        call handle_err(stat)
        stat = nf90_def_var(fid, "Y", nf90_float, yid, yvarid)
        call handle_err(stat)
        stat = nf90_def_var(fid, "time", nf90_float, rdimid, timeid)
        call handle_err(stat)
        stat = nf90_put_att(fid, timeid, 'units', 'seconds since 1987-03-26 00:00:00')
        call handle_err(stat)
        stat = nf90_enddef(fid)
        call handle_err(stat)
         
        stat = nf90_put_var(fid, xvarid, xaxis)
        call handle_err(stat)
        stat = nf90_put_var(fid, yvarid, yaxis)
        call handle_err(stat)
     endif

     varid=0
     
     do i=1, nvars
        if(trim(name)==trim(var(i)%nm)) then
           varid = var(i)%ID
           exit
        endif
     enddo
     
     if(varid==0) then
        if (nvars<50) then
           nvars = nvars+1
           var(nvars)%nm=name
           stat = nf90_redef(fid)
           call handle_err(stat)
           stat = nf90_def_var(fid, trim(name), nf90_float, (/xid,yid,rdimid/), var(nvars)%id)
           varid = var(nvars)%id
           call handle_err(stat)
           stat = nf90_enddef(fid)
           call handle_err(stat)
        else
           print *, 'Cannot write more than 50 vars'
           return
        endif
     endif
     
     stat = nf90_put_var(fid, varid, field, start=(/1,1,ntime1/))
     call handle_err(stat)

     stat = nf90_put_var(fid, timeid, (/ntime1/), start=(/ntime1/))
     call handle_err(stat)
    
     stat = nf90_sync(fid)
     call handle_err(stat)
    
   end subroutine write_nc_2d
   
  subroutine write_nc_2d_1(xaxis, yaxis, field, name, ntime)
     real, intent(in) :: xaxis(:), yaxis(:), field(:)
     character (len=*), intent(in) :: name
     integer, optional :: ntime
      
     !local vars
     integer :: xvarid, yvarid, stat, varid
     integer :: ntime1, i
     real :: tmp(size(xaxis),size(yaxis))
      
     ntime1 = 1
     if (present(ntime)) ntime1=ntime
      
     if (fid<0) then
        stat = nf90_create('output.nc', nf90_clobber, fid)
        call handle_err(stat)
        stat = nf90_def_dim(fid, "X", size(xaxis), xid)
        call handle_err(stat)
        stat = nf90_def_dim(fid, "Y", size(yaxis), yid)
        call handle_err(stat)
        stat = nf90_def_dim(fid, "time", nf90_unlimited, rdimid)
        call handle_err(stat)
         
        stat = nf90_def_var(fid, "X", nf90_float, xid, xvarid)
        call handle_err(stat)
        stat = nf90_def_var(fid, "Y", nf90_float, yid, yvarid)
        call handle_err(stat)
        stat = nf90_def_var(fid, "time", nf90_float, rdimid, timeid)
        call handle_err(stat)
        stat = nf90_put_att(fid, timeid, 'units', 'seconds since 1987-03-26 00:00:00')
        call handle_err(stat)
        stat = nf90_enddef(fid)
        call handle_err(stat)
         
        stat = nf90_put_var(fid, xvarid, xaxis)
        call handle_err(stat)
        stat = nf90_put_var(fid, yvarid, yaxis)
        call handle_err(stat)
     endif

     varid=0
     
     do i=1, nvars
        if(trim(name)==trim(var(i)%nm)) then
           varid = var(i)%ID
           exit
        endif
     enddo
     
     if(varid==0) then
        if (nvars<50) then
           nvars = nvars+1
           var(nvars)%nm=name
           stat = nf90_redef(fid)
           call handle_err(stat)
           stat = nf90_def_var(fid, trim(name), nf90_float, (/xid,yid,rdimid/), var(nvars)%id)
           varid = var(nvars)%id
           call handle_err(stat)
           stat = nf90_enddef(fid)
           call handle_err(stat)
        else
           print *, 'Cannot write more than 50 vars'
           return
        endif
     endif

     tmp = reshape(field,(/size(xaxis,1),size(yaxis,1)/))
     stat = nf90_put_var(fid, varid, tmp, start=(/1,1,ntime1/))
     call handle_err(stat)

     stat = nf90_put_var(fid, timeid, (/ntime1/), start=(/ntime1/))
     call handle_err(stat)
    
     stat = nf90_sync(fid)
     call handle_err(stat)
    
  end subroutine write_nc_2d_1

  subroutine handle_err(stat)
    integer, intent(in) :: stat
    if(stat /= nf90_noerr) then
      print *, trim(nf90_strerror(stat))
      stop "ERROR"
    end if
  end subroutine handle_err
  
end module netcdf_write
