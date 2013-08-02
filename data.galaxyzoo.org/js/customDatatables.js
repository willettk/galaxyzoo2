$(document).ready(function() {
	$('#example').dataTable( {					
		 "aaSorting": [[ 1, "asc" ]] // initially sorts column 2 (index starting at 0) ascending
	} );

} );