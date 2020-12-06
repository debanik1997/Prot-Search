import React, { useState, useEffect } from "react";
import "./ProteinEntry.css";
var parseString = require("xml2js").parseString;

export default function ProteinEntry(props) {
  const [isLoaded, setIsLoaded] = useState(false);
  const [protein, setProtein] = useState({});

  useEffect(() => {
    fetch("https://www.uniprot.org/uniprot/" + props.proteinID + ".xml")
      .then((response) => response.text())
      .then((xml) => {
        parseString(xml, function (err, result) {
          var entry = result.uniprot.entry[0];
          var proteinData = {
            proteinName: entry.name[0],
            protein: entry.protein[0].recommendedName[0].fullName[0],
            gene: entry.gene[0].name[0]._,
            organism: entry.organism[0].name[0]._,
            proteinID: props.proteinID,
          };
          setProtein(proteinData);
          setIsLoaded(true);
        });
      });
  }, [props.proteinID, props.offset]);

  if (!isLoaded) {
    return <></>;
  }
  return (
    <div className="proteinCard">
      <p>Name: {protein.proteinName}</p>
      <p>Protein: {protein.protein}</p>
      <p>Gene: {protein.gene}</p>
      <p>Organism: {protein.organism}</p>
      <p>
        More information at{" "}
        <a
          target="_blank"
          rel="noreferrer"
          href={"https://www.uniprot.org/uniprot/" + protein.proteinID}
        >
          Uniprot
        </a>
        .
      </p>
    </div>
  );
}
