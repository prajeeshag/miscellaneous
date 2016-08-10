module netcdf_write
  use netcdf
  
  implicit none
  private

  public :: write_nc

  interface write_nc
     module procedure write_nc_2d
  end interface write_nc
  
contains

  subroutine write_nc_2d(xaxis, yaxis, field, name)
    real, intent(in) :: xaxis(:), yaxis(:), field(:,:)
    character (len=*), intent(in) :: name
    !local vars
    integer :: fid, xid, yid, xvarid, yvarid, fieldid, stat
    
    stat = nf90_create('output.nc', nf90_clobber, fid)
    call handle_err(stat)
    stat = nf90_def_dim(fid, "X", size(xaxis), xid)
    call handle_err(stat)
    stat = nf90_def_dim(fid, "Y", size(yaxis), yid)
    call handle_err(stat)
    stat = nf90_def_var(fid, "X", nf90_float, xid, xvarid)
    call handle_err(stat)
    stat = nf90_def_var(fid, "Y", nf90_float, yid, yvarid)
    call handle_err(stat)
    stat = nf90_def_var(fid, trim(name), nf90_float, (/xid,yid/), fieldid)
    call handle_err(stat)
    stat = nf90_enddef(fid)
    call handle_err(stat)
     
    stat = nf90_put_var(fid, xvarid, xaxis)
    call handle_err(stat)
    stat = nf90_put_var(fid, yvarid, yaxis)
    call handle_err(stat)
    stat = nf90_put_var(fid, fieldid, field)
    call handle_err(stat)
    stat = nf90_close(fid)
    call handle_err(stat)
    
  end subroutine write_nc_2d

  subroutine handle_err(stat)
    integer, intent(in) :: stat
    if(stat /= nf90_noerr) then
      print *, trim(nf90_strerror(stat))
      stop "ERROR"
    end if
  end subroutine handle_err
  
end module netcdf_write
