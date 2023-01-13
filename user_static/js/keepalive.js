function session_keepalive() {
    //document.getElementById("logout-timer").innerHTML= "Attention: déconnection du planning dans " + logout_timer + " secondes si inactivité.";
    var xmlHttp = new XMLHttpRequest();
    xmlHttp.open( "GET", "./keepalive", false ); // false for synchronous request
    xmlHttp.send( null );
    setTimeout('session_keepalive()', 30000);
}

window.addEventListener('load', function(){
    setTimeout('session_keepalive()', 30000)
    document.getElementById("kal-warning").innerHTML = ""
})
