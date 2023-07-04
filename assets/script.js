var bgcColors = {
  NRPS: "",
  PKSI: "",
  "PKS-NRP_Hybrids": "",
  RiPPs: "",
  PKSother: "",
  Terpene: "",
  Others: "",
};

//background colors function
function actualizarColores(bgcKind) {
  var colorNRPS = document.getElementById("colorNRPS").value;
  var colorPKSI = document.getElementById("colorPKSI").value;
  var colorHybrids = document.getElementById("colorHybrids").value;
  var colorRiPPs = document.getElementById("colorRiPPs").value;
  var colorPKSother = document.getElementById("colorPKSothers").value;
  var colorTerpene = document.getElementById("colorTerpene").value;
  var colorOthers = document.getElementById("colorOthers").value;

  var bgcColors = {
    NRPS: colorNRPS,
    PKSI: colorPKSI,
    "PKS-NRP_Hybrids": colorHybrids,
    RiPPs: colorRiPPs,
    PKSother: colorPKSother,
    Terpene: colorTerpene,
    Others: colorOthers,
  };

  if (bgcKind != undefined) {
    // console.log(bgcKind);
    var pieChartImage = generarPieChart(bgcKind, bgcColors);
    // nodo.style("background-image", "url(" + pieChartImage + ")");
    return pieChartImage;
    // Establecer el estilo de Cytoscape.js para el nodo actual
  } else {
    // console.log(bgcColors);
    return bgcColors;
  }
}

// // Actualizar los colores por defecto al cargar la página
// window.addEventListener("DOMContentLoaded", function () {
//   actualizarColores();
//   // handleFile();
// });

function getElementsInsideNode(nodeId, data, kind) {
  var currentData = data.find(function (diccionario) {
    return diccionario.GCF === nodeId;
  });
  if (currentData != undefined) {
    var elementosBGC = [];
    var estructuraCorregida = currentData.GCF_data.replace(/'/g, '"').replace(
      /nan/g,
      "null"
    );
    let diccionarioBGCs = JSON.parse(estructuraCorregida);
    var elementosBGC = diccionarioBGCs.map(function (diccionario) {
      return diccionario[kind];
    });

    return elementosBGC;
  }
}

// Función para generar el piechart según la data del nodo
function generarPieChart(data, colores) {
  // console.log(data);
  var canvas = document.createElement("canvas");
  var ctx = canvas.getContext("2d");
  var total = data.reduce(function (sum, item) {
    return sum + item.value;
  }, 0);

  var startAngle = 0;
  data.forEach(function (item, index) {
    var sliceAngle = (item.value / total) * 2 * Math.PI;

    ctx.beginPath();
    ctx.moveTo(canvas.width / 2, canvas.height / 2);
    ctx.arc(
      canvas.width / 2,
      canvas.height / 2,
      canvas.width / 2,
      startAngle,
      startAngle + sliceAngle
    );
    ctx.fillStyle = colores[item.BGCKind];
    ctx.fill();

    startAngle += sliceAngle;
  });

  var pieChartImage = canvas.toDataURL(); // Obtener imagen como base64
  return pieChartImage;
}

function handleFile() {
  var networkFile = document.getElementById("fileInput");
  const file = networkFile.files[0];

  const reader = new FileReader();
  reader.onload = function (event) {
    const fileContent = event.target.result;
    var jsonData = JSON.parse(fileContent);

    makeNetwork(jsonData);
  };
  // console.log(jsonData);
  reader.readAsText(file);

  makeNetwork = (data1) => {
    var elements = {
      nodes: [],
      edges: [],
    };

    // Construir los nodos
    data1.forEach((data) => {
      elements.nodes.push({
        data: {
          id: data.MFs,
          style: { "background-color": "#80CB80", shape: "triangle" },
        },
      });
    });

    data1.forEach((data) => {
      elements.nodes.push({
        data: {
          id: data.GCF,
          style: {
            "background-color": "#DDAE70",
          },
        },
      });
    });

    // Construir las aristas
    data1.forEach((data) => {
      var source = data.MFs;
      var target = data.GCF;
      var weight = data.jaccard_score;
      var weightFormatted = weight.toFixed(1); // Limitar a 1 dígito después de la coma

      elements.edges.push({
        data: { source: source, target: target, weight: weightFormatted },
      });
    });
    // Construir los nodos
    data1.forEach((data) => {
      elements.nodes.push({
        data: {
          id: data.MFs,
        },
      });
    });

    // Crear la visualización de la red
    var cy = cytoscape({
      container: document.getElementById("cy"),
      maxZoom: 3,
      minZoom: 0.5,
      elements: elements,
      // console.log(elements),
      style: [
        {
          selector: "node",
          // group: "nodes",
          style: {
            "background-color": "grey",

            // "url(Ejemplo.svg)",
            "border-width": 1.5,
            "border-style": "solid",

            shape: function (node) {
              var nodeName = node.data("id");
              if (nodeName.includes("GCF")) {
                return "ellipse"; // Color para nodos GCF
              } else if (nodeName.includes("MF")) {
                return "triangle"; // Color para nodos MF
              } else {
                return "ellipse";
              }
            },
            label: "data(id)",
            // width: 20, // Ajustar el tamaño del nodo
            "font-size": 10, // Ajustar el tamaño de la fuente
          },
        },

        {
          selector: "edge",
          style: {
            "curve-style": "straight",
            "line-color": "#ccc",
            label: "data(weight)",
            "text-valign": "center",
            "text-halign": "center",
            "font-size": 12,
            "source-endpoint": "outside-to-node",
            "target-endpoint": "outside-to-node",
          },
        },
      ],
      layout: {
        name: "cose",
        idealEdgeLength: 100, // Ajustar la longitud ideal de las aristas
        nodeOverlap: 200, // Ajustar el valor de superposición de nodos
        padding: 100, // Ajustar el espacio de relleno alrededor de la red
        randomize: false, // Desactivar la aleatorización del layout
        componentSpacing: 100, // Ajustar el espacio entre componentes
        //   rows: 10,
      },
    });

    // trabajar esto!!
    // searchInput.addEventListener("input", function (event) {
    //   var searchText = event.target.value.toLowerCase();

    //   cy.nodes().forEach(function (node) {
    //     var nodeId = node.id().toLowerCase();

    //     if (nodeId.includes(searchText)) {
    //       console.log(nodeId);
    //       node.addClass("highlight");
    //       console.log(node);
    //     } else {
    //       node.removeClass("highlight");
    //     }
    //   });
    // });

    // Obtén todos los nodos del grafo
    var nodos = cy.nodes();

    // Itera sobre cada nodo y aplica el estilo del piechart individual
    nodos.forEach(function (nodo) {
      var nodeId = nodo.id();
      if (nodeId.includes("GCF")) {
        var data = getElementsInsideNode(nodeId, data1, "class"); // Obtiene la data correspondiente al nodo actual
        // console.log(data);

        // Contar los elementos y almacenar su frecuencia en un objeto
        const elementCount = {};
        data.forEach((element) => {
          elementCount[element] = (elementCount[element] || 0) + 1;
        });

        // Ordenar las claves del objeto
        const sortedKeys = Object.keys(elementCount).sort();

        // Crear el arreglo final en formato de diccionario
        const bgcKind = sortedKeys.map((key) => ({
          BGCKind: key,
          value: elementCount[key],
        }));

        nodo.style({
          // "background-color": "grey",
          "background-image": actualizarColores(bgcKind),
        });
      }
    });

    // Mostrar elementos contenidos en el nodo al hacer hover
    var popover = null; // Variable para almacenar el popover actual
    cy.on("mouseover", "node", function (event) {
      var node = event.target;
      var nodeId = node.id();
      var elementsInsideNode = getElementsInsideNode(nodeId, data1, "BGC"); // Función para obtener los elementos contenidos en el nodo

      // Crear el contenido del popover
      var popoverContent = "<p>ID del nodo: " + nodeId + "</p>";
      elementsInsideNode.forEach((element) => {
        popoverContent += "<p>" + element + "</p>";
      });

      // Eliminar el popover existente si hay uno
      if (popover !== null) {
        popover.remove();
        popover = null;
      }

      // Crear el popover
      popover = document.createElement("div");
      popover.classList.add("popover");
      popover.classList.add("text-center");
      popover.innerHTML = popoverContent;

      // Posicionar el popover sobre el nodo
      var nodePosition = node.renderedPosition();
      popover.style.left = nodePosition.x + "px";
      popover.style.top = nodePosition.y + "px";

      // Agregar el popover al contenedor de la visualización
      var cyContainer = document.getElementById("cy");
      cyContainer.appendChild(popover);
    });

    cy.on("mouseout", "node", function (event) {
      // Eliminar el popover cuando se sale del nodo
      if (popover !== null) {
        popover.remove();
        popover = null;
      }
    });

    cy.on("mousedown", function (event) {
      // Ocultar el popover al hacer clic en cualquier lugar fuera de los nodos
      if (popover !== null) {
        popover.remove();
        popover = null;
      }
    });

    function getRandomOffset(index, total) {
      // Calcular el ángulo basado en el índice y el total de nodos
      var angle = (2 * Math.PI * index) / total;

      // Definir el radio del círculo
      var radius = 50; // Ajusta el valor del radio según tus necesidades

      // Calcular las coordenadas x e y en función del ángulo y el radio
      var x = radius * Math.cos(angle);
      var y = radius * Math.sin(angle);

      return { x: x, y: y };
    }

    cy.on("tap", "node", function (event) {
      var gcfNode = event.target;
      var currentStyle = gcfNode.style();
      var bgColor = currentStyle["background-color"];
      var bgcs_edges = obtenerBGCsDeGCF(gcfNode.id());

      // var bgc_info;

      var bgcs_class = getElementsInsideNode(gcfNode.id(), data1, "class");
      var bgcs_name = getElementsInsideNode(gcfNode.id(), data1, "BGC");
      var longitud = Math.min(bgcs_class.length, bgcs_name.length);

      // Crear un nuevo array para almacenar los resultados combinados
      var bgc_info = [];

      // Combinar los elementos de ambos arrays
      for (var i = 0; i < longitud; i++) {
        bgc_info.push({
          BGC: bgcs_name[i],
          class: bgcs_class[i],
        });
      }

      // Imprimir el resultado combinado
      // console.log(resultado);

      if (bgc_info) {
        if (gcfNode.style("width") != "200px") {
          // Cambiar el tamaño del nodo al hacer clic
          gcfNode.style("width", "200px");
          gcfNode.style("height", "200px");
          gcfNode.style("background-color", "red");
          gcfNode.style("background-image-opacity", 0);
          gcfNode.style("background-opacity", 0.1);
          // gcfNode.panna;
          gcfNode.lock();

          bgc_info.forEach(function (bgc, index) {
            // Objeto auxiliar para realizar un seguimiento de los nodos BGC agregados
            var addedNodes = {};

            // Agregar los nodos BGCs con la posición deseada
            // var randomOffset = getRandomOffset();
            var randomOffset = getRandomOffset(index, bgc_info.length);

            // console.log(randomOffset);
            // console.log(gcfNode.position());
            var newPosition = {
              x: gcfNode.position().x + randomOffset.x,
              y: gcfNode.position().y + randomOffset.y,
            };

            // console.log(newPosition);
            // var uniqueId = bgc.source + "-" + bgc.target; // Generar un ID único para el nodo BGC
            if (!addedNodes[bgc.BGC]) {
              var bgcColors = actualizarColores();
              // console.log(bgcColors);
              // var bgcClass = bgc.class;
              // console.log(bgc);
              let nodeId = bgc.BGC + "/" + gcfNode.id();
              backgroundColor = bgcColors[bgc.class];
              // console.log(backgroundColor);
              // console.log(backgroundColor);
              try {
                cy.add({
                  parent: gcfNode.id(),
                  data: { id: nodeId },
                  // position: { parent: gcfNode.position() },
                  position: {
                    x: newPosition.x,
                    y: newPosition.y,
                    //   // newPosition,
                  },
                  style: {
                    "background-color": backgroundColor, // Ajusta el valor del color según tus necesidades
                    label: bgc.BGC,
                  },
                  // classes: "bgc-node",
                });
              } catch {
                console.log("node already added");
              }
            }
          });

          bgcs_edges.forEach(function (bgc) {
            let source = bgc.source + "/" + gcfNode.id();
            let target = bgc.target + "/" + gcfNode.id();
            // let weight = bgc.weight[0:2]
            try {
              cy.add({
                data: {
                  source: source,
                  target: target,
                  weight: bgc.weight.toString().substring(0, 3),
                },
                classes: "bgc-edge",
              });
            } catch {
              console.log("cannot create something");
            }
          });
        } else {
          gcfNode.style("height", "30px");
          gcfNode.style("width", "30px");
          gcfNode.style("background-image-opacity", 1);
          gcfNode.style("background-color", bgColor);
          gcfNode.unlock();

          bgcs_edges.forEach(function (bgc, index) {
            let source = bgc.source + "/" + gcfNode.id();
            let target = bgc.target + "/" + gcfNode.id();
            var nodeSource = cy.getElementById(source);
            var nodeTarget = cy.getElementById(target);
            nodeSource.remove();
            nodeTarget.remove();
          });
        }
      }
      // Definir el rango de movimiento relativo
    });

    // Función para obtener los BGCs asociados a un nodo de GCF específico
    function obtenerBGCsDeGCF(gcfNode) {
      // Aquí puedes realizar la lógica para obtener los datos de los BGCs asociados al GCF específico
      // Puedes usar gcfNodeId para filtrar los datos necesarios

      if (gcfNode.includes("GCF")) {
        // Ejemplo de datos de BGCs
        try {
          // console.log(data1);
          // console.log(networkData);
          const elementoEncontrado = data1.find(
            (elemento) => elemento.GCF === gcfNode
          );

          if (elementoEncontrado) {
            const gcfEdges = elementoEncontrado.GCF_edges;
            // console.log(gcfEdges[0]);

            const resultados = [];

            for (let elemento of gcfEdges) {
              const [source, target, weight] = elemento
                .substring(1, elemento.length - 1)
                .split(",");
              resultados.push({
                source: source.trim(),
                target: target.trim(),
                weight: parseFloat(weight.trim()),
              });
            }

            // console.log(resultados);
            return resultados;
          } else {
            console.log("Element not found");
          }

          // var bgcsData = data1[gcfNode].BGCs.map(function (bgc) {
          //   return {
          //     source: bgc.source,
          //     target: bgc.target,
          //     weight: bgc.weight,
          //   };
          // });
          return bgcsData;
        } catch {
          console.log("This node has no data");
        }
      }
      // Retorna los datos de los BGCs como un arreglo
    }

    var moveRange = 70; // Valor máximo de movimiento en cualquier dirección

    // Escuchar el evento 'dragfree' de los nodos
    // bgcs.forEach(function (bgc, index) {
    // var node = cy.getElementById(bgc.source);
    cy.on("drag", "node", function (event) {
      var node = event.target;
      var nodePosition = node.position();
      var parentNode = node.id().split("/")[1];
      // console.log(parentNode);
      // var nodeID = node.id();
      // console.log(nodeID);

      // var currentGCF = buscarGCF(node.id(), datosRed);
      // console.log(currentGCF);
      var GCFnode = cy.getElementById(parentNode); // Obtener el nodo por su ID
      var GCFposition = GCFnode.position(); // Obtener la posición del nodo
      // var GCFposition = obtenerBGCsDeGCF(currentGCF.position());
      // console.log(GCFposition);

      if (node.id().includes("MF")) {
        console.log("MF node");
      } else if (node.id().includes("GCF") && !node.id().includes("/")) {
        console.log("GCF node");
      } else {
        console.log("bgc node");
        // console.log(nodePosition);

        // Calcular los límites teniendo en cuenta el rango de movimiento relativo
        var minX = GCFposition.x - moveRange;
        var maxX = GCFposition.x + moveRange;
        var minY = GCFposition.y - moveRange;
        var maxY = GCFposition.y + moveRange;

        // Verificar si la posición del nodo está fuera de los límites
        var newPosX = Math.max(minX, Math.min(maxX, nodePosition.x));
        var newPosY = Math.max(minY, Math.min(maxY, nodePosition.y));

        // Establecer la nueva posición del nodo dentro de los límites
        node.position({ x: newPosX, y: newPosY });
      }
    });
    // });

    // arreglar esto y estamos operativos!!!
    function buscarGCF(string, red) {
      for (var gcf in red) {
        var bgcs = red[gcf].BGCs;

        for (var i = 0; i < bgcs.length; i++) {
          if (bgcs[i].source === string || bgcs[i].target === string) {
            return gcf;
          }
        }
      }

      return null; // Retorna null si el string no se encuentra en ninguna GCF
    }
  };
}

// onload(handleFile());
