self: super: {
  # Override netcdf build setttings for m1 osx
	netcdf = super.netcdf.overrideAttrs (
    old: {
	    configureFlags = [
	         "--enable-netcdf-4"
	         "--enable-dap"
	         "--enable-shared"
	         "--disable-dap-remote-tests"
	         "--disable-unit-tests"
	         "--disable-testsets"
	         "--disable-examples"
	         "--disable-filter-testing"
	       ];
    }
	);
}
