import React, { useState } from "react";
import { Redirect } from "react-router-dom";
import LoadingDNA from "./LoadingDNA";
import Checkbox from "./Checkbox";
import "./SearchPage.css";

export default function SearchPage() {
  const [searchText, setSearchText] = useState("");
  const [maxMismatches] = useState(4);
  const [gapsAllowed, setGapsAllowed] = useState(false);
  const [results, setResults] = useState([]);
  const [redirect, setRedirect] = useState(false);

  const handleSearch = () => {
    const form = {
      max_mismatches: maxMismatches,
      pattern: searchText,
      gaps_allowed: gapsAllowed,
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
              setRedirect(true);
            }
          });
        } else {
          alert("Error with search");
        }
      })
      .catch(() => {
        alert("Error with search");
      });
  };

  if (redirect) {
    return (
      <Redirect
        push
        to={{
          pathname: "/results",
          results: results,
        }}
      />
    );
  }

  return (
    <>
      <div className="searchPage">
        <h1>Welcome to ProtSearch</h1>
        <LoadingDNA />
        <input
          name="text"
          type="text"
          placeholder="Type a DNA read here..."
          value={searchText}
          onChange={(e) => setSearchText(e.target.value)}
        />
        <button className="searchButton" onClick={handleSearch}>
          Search
        </button>
        <div>
          <label>
            <Checkbox
              checked={gapsAllowed}
              onChange={(e) => setGapsAllowed(e.target.checked)}
            />
            <span className="checkboxText" style={{marginLeft: 10}}>Insertions and Deletions allowed?</span>
          </label>
        </div>
      </div>
    </>
  );
}
