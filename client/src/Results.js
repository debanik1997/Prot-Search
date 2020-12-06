import React, { useState, useEffect } from "react";
import ProteinEntry from "./ProteinEntry";

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
    return <h1 className="results">Searching...</h1>;
  }

  return (
    <>
      <h1 className="results">Results</h1>
      {results.map((result) => {
        return (
          <ProteinEntry
            key={result[0]}
            proteinID={result[0]}
            offset={result[1]}
          />
        );
      })}
    </>
  );
}
