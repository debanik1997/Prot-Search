import React, { useState } from "react";
import { Redirect } from "react-router-dom";
import "./SearchPage.css"

export default function SearchPage() {
  const [searchText, setSearchText] = useState("");
  const [maxMismatches] = useState(4);
  const [gapsAllowed] = useState(false);
  const [redirect, setRedirect] = useState(false);

  const handleSearch = () => {
    setRedirect(true);
  };

  if (redirect) {
    return (
      <Redirect
        push
        to={{
          pathname: "/results",
          searchText: searchText,
          maxMismatches: maxMismatches,
          gapsAllowed: gapsAllowed,
        }}
      />
    );
  }

  return (
    <div>
      <h1>Welcome to ProtSearch</h1>
      <input
        name="text"
        type="text"
        placeholder="Search"
        value={searchText}
        onChange={(e) => setSearchText(e.target.value)}
      />
      <button onClick={handleSearch}>Search</button>
    </div>
  );
}
