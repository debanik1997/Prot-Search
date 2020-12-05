import React, { useState, useEffect } from "react";

export default function Results(props) {
  const [isLoaded, setIsLoaded] = useState(false);
  const [results, setResults] = useState([]);

  useEffect(() => {
    if (props.location && props.location.results) {
      setIsLoaded(true);
      setResults(props.location.results);
    }
  }, [props.location]);

  if (!isLoaded) {
    return<h1 className="results">Searching...</h1>
  }

  return results.map((result) => {
    return (
      <>
        <h1 className="results">Results</h1>
        <div
          key={result[0]}
          style={{
            height: "120px",
            width: "200px",
            border: "1px solid #CECECE",
            boxSizing: "border-box",
            borderRadius: "0px 0px 6px 6px",
            borderTop: "4px solid #F54CEE",
            marginBottom: 16,
            marginRight: 10,
            marginLeft: 10,
            cursor: "auto",
            background: "#FFFFFF",
          }}
        >
          <p>Protein ID: {result[0]}</p>
          <p>Offset: {result[1]}</p>
          <a
            target="_blank"
            rel="noreferrer"
            href={"https://www.uniprot.org/uniprot/" + result[0]}
          >
            More information
          </a>
        </div>
      </>
    );
  });
}
