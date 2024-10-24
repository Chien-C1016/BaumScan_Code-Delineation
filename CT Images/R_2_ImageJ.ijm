// R 2 ImageJ is for converting the semi-automatic tree rings
// into interactive ImageJ platform
// In such way, it is much easier to check which tree ring is wrongly detected.

// Ask for the CSV file
path = File.openDialog("Choose a CSV File");
content = File.openAsString(path);
lines = split(content, "\n");

// Initialize arrays to store labels and their corresponding coordinates as strings
labelMap = newArray();
coordsMap = newArray();

// Process each line, starting from 1 to skip the header
for (i = 1; i < lines.length; i++) {
    values = split(lines[i], ",");
    if (values.length == 3) { // Ensure there are enough values per line
        x = values[0];
        y = values[1];
        label = values[2];
        
        // Manual search for label index
        index = -1;
        for (j = 0; j < lengthOf(labelMap); j++) {
            if (labelMap[j] == label) {
                index = j;
                break;
            }
        }
        
        if (index == -1) {
            // New label, add it to the array
            labelMap = Array.concat(labelMap, label);
            coordsMap = Array.concat(coordsMap, x + "," + y); // Initialize with first coordinate pair
        } else {
            // Existing label, append coordinates
            coordsMap[index] = coordsMap[index] + "|" + x + "," + y;
        }
    }
}

// Prepare the ROI Manager
roiManager("Reset");

// Create and add polygon ROIs to ROI Manager for each label
for (i = 0; i < lengthOf(labelMap); i++) {
    coordPairs = split(coordsMap[i], "|");
    xCoords = newArray();
    yCoords = newArray();
    for (j = 0; j < lengthOf(coordPairs); j++) {
        coords = split(coordPairs[j], ",");
        xCoords = Array.concat(xCoords, parseFloat(coords[0]));
        yCoords = Array.concat(yCoords, parseFloat(coords[1]));
    }
    makeSelection("polygon", xCoords, yCoords);
    roiManager("Add");
    roiManager("select", roiManager("Count")-1);
    roiManager("Rename", "Ring-" + labelMap[i]);
}

print("Done processing CSV file.");

