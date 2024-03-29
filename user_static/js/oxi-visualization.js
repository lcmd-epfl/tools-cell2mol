

function toggleStrVisInteraction(enableStrInteraction) {
  if (enableStrInteraction) {
    // enable interaction here
    $("#str-overlay")
      .css("display", "none")
      .css("-webkit-touch-callout", "auto")
      .css("-webkit-user-select", "auto")
      .css("-khtml-user-select", "auto")
      .css("-moz-user-select", "auto")
      .css("-ms-user-select", "auto")
      .css("user-select", "auto");
  } else {
    // disable interaction here
    $("#str-overlay")
      .css("display", "table")
      .css("-webkit-touch-callout", "none")
      .css("-webkit-user-select", "none")
      .css("-khtml-user-select", "none")
      .css("-moz-user-select", "none")
      .css("-ms-user-select", "none")
      .css("user-select", "none");
  }
}

//function jsmolCrystal(data, parentHtmlId, appletName, supercellOptions) {
function jsmolCrystal(data, ucell, parentHtmlId, appletName, supercellOptions, atmCon) {
  var parentDiv = document.getElementById(parentHtmlId);
  var the_width = parentDiv.offsetWidth - 5;
  var the_height = parentDiv.offsetHeight - 5;

  var Info = {
    width: the_width,
    height: the_height,
    debug: false,
    color: "#FFFFFF",
    use: "HTML5",
    j2sPath: "../user_static/js/jsmol/j2s",
    serverURL: "../user_static/js/jsmol/php/jsmol.php",
    console: "jmolApplet_infodiv"
  };

  var jsmolStructureviewer = Jmol.getApplet(appletName, Info);

  if (supercellOptions === undefined) {
    var loadingScript =
      'color cpk; load INLINE "' + data + '"'+"{ijk i'j'k' 0}"+' centroid unitcell "' + ucell + '"; unitcell true; select all; hideNotSelected = true; zoom 50;'
      + atmCon + ' select all ;';
  } else {
    var loadingScript =
      'color cpk; load INLINE "' +
      data +
      '" {' +
      supercellOptions[0] +
      " " +
      supercellOptions[1] +
      " " +
      supercellOptions[2] +
      "} centroid unitcell \""
      ucell + '"';
  }

  // set unit cell data
  //loadingScript += "; unitcell \""+ucell+"\"";
  loadingScript += "; frame all; hide none";

  //Jmol.script(jsmolStructureviewer, loadingScript)
    
  //draw x,y,z axes instead of a,b,c vectors as default
  loadingScript +=
    '; axes off; draw xaxis ">X" vector {0 0 0} {2 0 0} color red width 0.15; draw yaxis ">Y" vector {0 0 0} {0 2 0} color green width 0.15; draw zaxis ">Z" vector {0 0 0} {0 0 2} color blue width 0.15';

  loadingScript += "; wireframe 0.1; spacefill 23%";
  loadingScript += "; unitcell primitive";

  //Sets the unit cell line diameter in Angstroms
  loadingScript += "; unitcell 2";

  // antialiasDisplay ON
  loadingScript += "; set antialiasDisplay on";

  //Zooms to the setting that fills the screen with the currently displayed atoms
  loadingScript += "; set zoomLarge false";

  //remove JSmol logo
  loadingScript += "; set frank off";

  //setup labels
  loadingScript += "; set labeloffset 2 2";
  loadingScript += "; set fontSize 16";

  loadingScript += "select all"


  Jmol.script(jsmolStructureviewer, loadingScript);

  //parentDiv.innerHTML = Jmol.getAppletHtml(jsmolStructureviewer);
  //$("#" + parentHtmlId).html(Jmol.getAppletHtml(jsmolStructureviewer));

  return jsmolStructureviewer;
}

var cellLine = "; unitcell 2";

function toggleRotation(viewer) {
  if ($("#spin-input").is(":checked")) {
    var jmolscript = "spin on";
  } else {
    var jmolscript = "spin off";
  }
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

function cmpVisibilityUpdate(viewer) {                                                                                                
  let visible_cmps = [];                                                                                                        
  let visible_cmps_str = "";                                                                                                    
  let x = document.getElementById("structure-chooser");                                                                         
  let i;                                                                                                                        
  visible_cmps_str += "select   ";
  for (i = 0; i < x.length ;i++) {                                                                                              
    if (x.elements[i].checked) {                                                                                                
      visible_cmps.push(i);                                                                                                         
      visible_cmps_str += (" " + x.elements[i].value + " or");                                                                                             
    }                                                                                                                           
  }                                                                                                                             
  let jmolscript = visible_cmps_str.slice(0,-2);
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}                                                

function showUnpacked(viewer) {
  var jmol_list_pos = document.getElementById("atm_pos").value.split('$')[0];
  if ($("#unpacked-input").is(":checked")) {
    //var jmolscript = jmol_list_pos+ "; unitcell false";
    var jmolscript = "model " + jmol_list_pos + "; unitcell false";
  } else {
    //var jmolscript = "select all; unitcell true";
    var jmolscript = "model 0; unitcell true";
  }
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

function c2mButton(viewer) {
  var jmol_list_pos = document.getElementById("atm_pos").value.split('$')[0];
  if ($("#unpacked-input").is(":checked")) {
    document.getElementById("atm_pos").style.display="block";  
    document.getElementById("label_pos").style.display="block";  
    document.getElementById("downloadBtn").style.display="block";  
    var jmolscript = "model " + jmol_list_pos + "; unitcell false";
  } else {
    document.getElementById("atm_pos").style.display="none";
    document.getElementById("label_pos").style.display="none";
    document.getElementById("downloadBtn").style.display="none";
    //var jmolscript = "select all; unitcell true";
    var jmolscript = "model 0; unitcell true";
  }
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

function showBonds(viewer) {
  if ($("#bonds-input").is(":checked")) {
    var jmolscript = "wireframe 0.1";
  } else {
    var jmolscript = "wireframe off";
  }
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

function showPacked(viewer) {
  var nx = $("#nx").val();
  var ny = $("#ny").val();
  var nz = $("#nz").val();

  if ($("#packed-input").is(":checked")) {
    var jmolscript =
      "save orientation 0; load '' {" +
      nx +
      " " +
      ny +
      " " +
      nz +
      "} packed; unitcell primitive; restore orientation 0" +
      jsmolDrawAxes(viewer) +
      cellLine +
      "; " +
      showLabels(viewer) +
      "; " +
      showBonds(viewer);
  } else {
    var jmolscript =
      "save orientation 0; load '' {" +
      nx +
      " " +
      ny +
      " " +
      nz +
      "}; unitcell primitive; restore orientation 0" +
      jsmolDrawAxes(viewer) +
      cellLine +
      "; " +
      showLabels(viewer) +
      "; " +
      showBonds(viewer);
  }
  $("#spin-input").prop("checked", false);
  $("#spheres-input").prop("checked", false);
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

function showLabels(viewer) {
  if ($("#labels-input").is(":checked")) {
    var jmolscript = "label %e";
  } else {
    var jmolscript = "label off";
  }
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

function labelOxStates(viewer, atom_indices, labeltext) {
  var jmolscript = "";
  var counter = 0;
  const elements = labeltext.length;

  if ($("#labels-input").is(":checked")) {
    for (const [index, label] of labeltext.entries()) {
      if (counter < elements) {
        var num = atom_indices[index] + 1;
        jmolscript += "select atomno=" + num + "; ";
        jmolscript += "label " + '"' + label + '"' + "; ";
        counter++;
      } else {
        var num = atom_indices[index] + 1;
        jmolscript += "select atomno=" + num + "; ";
        jmolscript += "label " + '"' + label + '" select none';
      }
    }
  } else {
    for (const [index, label] of labeltext.entries()) {
      if (counter < elements) {
        var num = atom_indices[index] + 1;
        jmolscript += "select atomno=" + num + "; ";
        jmolscript += "label off; ";
        counter++;
      } else {
        var num = atom_indices[index] + 1;
        jmolscript += "select atomno=" + num + "; ";
        jmolscript += "label off; select none";
      }
    }
  }
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

function showSpheres(viewer) {
  if ($("#spheres-input").is(":checked")) {
    var jmolscript = "spacefill on";
  } else {
    var jmolscript = "spacefill 23%";
  }
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

function jsmolDrawAxes(viewer) {
  var e = document.getElementById("axesMenu");
  var selectedAxes = e.options[e.selectedIndex].value;
  switch (selectedAxes) {
    case "xyz":
      var jmolscript =
        "; axes off; draw xaxis '>X' vector {0 0 0} {2 0 0} color red width 0.15; draw yaxis '>Y' vector {0 0 0} {0 2 0} color green width 0.15; draw zaxis '>Z' vector {0 0 0} {0 0 2} color blue width 0.15";
      break;
    case "abc":
      var jmolscript =
        "; draw xaxis delete; draw yaxis delete; draw zaxis delete; set axesMode 2; axes 5";
      break;
    case "noaxes":
      var jmolscript =
        "; draw xaxis delete; draw yaxis delete; draw zaxis delete; axes off";
  }
  Jmol.script(eval(viewer), jmolscript);
  return jmolscript;
}

// function jsmolSupercell(viewer) {
//   var nx = $("#nx").val();
//   var ny = $("#ny").val();
//   var nz = $("#nz").val();
//   $("#spin-input").prop("checked", false);
//   $("#spheres-input").prop("checked", false);
//   $("#packed-input").prop("checked", false);
//   var jmolscript =
//     "save orientation 0; load '' {" +
//     nx +
//     " " +
//     ny +
//     " " +
//     nz +
//     "}; unitcell primitive; restore orientation 0" +
//     jsmolDrawAxes(viewer) +
//     cellLine +
//     "; " +
//     showLabels(viewer) +
//     "; " +
//     showBonds(viewer);
//   Jmol.script(eval(viewer), jmolscript);
// }

// function jsmol222cell(viewer) {
//   $("#spin-input").prop("checked", false);
//   $("#spheres-input").prop("checked", false);
//   $("#packed-input").prop("checked", false);
//   // reset nx, ny, nz to 2,2,2
//   $("#nx").val(2);
//   $("#ny").val(2);
//   $("#nz").val(2);
//   Jmol.script(
//     eval(viewer),
//     "save orientation 0; load '' {2 2 2}; unitcell primitive; restore orientation 0" +
//       jsmolDrawAxes(viewer) +
//       cellLine +
//       "; " +
//       showLabels(viewer) +
//       "; " +
//       showBonds(viewer)
//   );
// }

function centerXaxis(viewer) {
  Jmol.script(eval(viewer), "moveto 1 axis x");
}

function centerYaxis(viewer) {
  Jmol.script(eval(viewer), "moveto 1 axis y");
}

function centerZaxis(viewer) {
  Jmol.script(eval(viewer), "moveto 1 axis z");
}

function viewerPostLoad(viewer) {
    Jmol.script(eval(viewer), "set displayCellParameters FALSE; set antialiasDisplay on");
}
function showCompounds(viewer, visible) {
    Jmol.script(eval(viewer), "frame [" + visible + "]");
}

$.fn.bindFirst = function(name, fn) {
  var elem, handlers, i, _len;
  this.bind(name, fn);
  for (i = 0, _len = this.length; i < _len; i++) {
    elem = this[i];
    handlers = jQuery._data(elem).events[name.split(".")[0]];
    handlers.unshift(handlers.pop());
  }
};

function enableDoubleTap(element, callback, ignoreOnMove) {
  /* Enable double-tap event for phones */
  element.dbltapTimeout = undefined;
  element.shortTap = false;

  var preventOnMove = false;
  if (typeof ignoreOnMove !== "undefined") {
    preventOnMove = ignoreOnMove;
  }

  // Manual detect of double tap
  //element.addEventListener('touchend', function(event) {
  $(element).bindFirst("touchend", function(event) {
    if (typeof element.dbltapTimeout !== "undefined") {
      // start disabling any timeout that would reset shortTap to false
      clearTimeout(element.dbltapTimeout);
    }
    if (element.shortTap) {
      // if here, there's been another tap a few ms before
      // reset the variable and do the custom action
      element.shortTap = false;
      event.preventDefault();
      event.stopImmediatePropagation();
      callback();
    } else {
      if (event.targetTouches.length != 0) {
        // activate this only when there is only a finger
        // if more than one finger is detected, cancel detection
        // of double tap
        if (typeof element.dbltapTimeout !== "undefined") {
          // disable the timeout
          clearTimeout(element.dbltapTimeout);
          element.shortTap = false;
        }
        return;
      }
      // If we are here, no tap was recently detected
      // mark that a tap just happened, and start a timeout
      // to reset this
      element.shortTap = true;

      element.dbltapTimeout = setTimeout(function() {
        // after 500ms, reset shortTap to false
        element.shortTap = false;
      }, 500);
    }
  });
  element.addEventListener("touchcancel", function(event) {
    if (typeof element.dbltapTimeout !== "undefined") {
      // disable the timeout if the touch was canceled
      clearTimeout(element.dbltapTimeout);
      element.shortTap = false;
    }
  });
  if (!preventOnMove) {
    element.addEventListener("touchmove", function(event) {
      if (typeof element.dbltapTimeout !== "undefined") {
        // disable the timeout if the finger is being moved
        clearTimeout(element.dbltapTimeout);
        element.shortTap = false;
      }
    });
  }
}




