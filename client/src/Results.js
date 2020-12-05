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
    return <div>Searching...</div>;
  }

  return results.map((result) => {
    return (
      <div key={result[0]}>
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
    );
  });
}
