import React, { useState } from "react";
export default function SearchPage() {
  const [searchText, setSearchText] = useState("");

  const handleSearch = () => {
    console.log(searchText);
  };

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
