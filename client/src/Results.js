import React, { useState, useEffect } from "react";

export default function Results(props) {
  const [isLoaded, setIsLoaded] = useState(false);
  const [results, setResults] = useState([]);

  useEffect(() => {
    if (props.location && props.location.searchText) {
      const form = {
        max_mismatches: props.location.maxMismatches,
        pattern: props.location.searchText,
        gaps_allowed: props.location.gapsAllowed,
      };
      fetch("/api/protein", {
        method: "post",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(form),
      })
        .then((response) => {
          if (response.ok) {
            response.json().then((data) => {
              if (data.match) {
                setResults(data.match);
                setIsLoaded(true);
              }
            });
          } else {
            alert("Error with search");
          }
        })
        .catch(() => {
          alert("Error with search");
        });
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
