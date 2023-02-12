$(document).on("click", ".go-map", function(e) {
  e.preventDefault();
  $el = $(this);
  var lat = $el.data("lat");
  var long = $el.data("long");
  var pin = $el.data("pin");
  $($("#nav a")[0]).tab("show");
  Shiny.onInputChange("goto", {
    lat: lat,
    lng: long,
    pin: pin,
    nonce: Math.random()
  });
});