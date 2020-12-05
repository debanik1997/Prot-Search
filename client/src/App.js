import React, { useState } from "react";
import SearchBar from "material-ui-search-bar";
import "./App.css";

function App() {
  const [searchText, setSearchText] = useState('');

  return (
    <div className="content">
      <SearchBar
        value={searchText}
        onChange={(newSearchText) => setSearchText(newSearchText)}
        onRequestSearch={() => console.log(searchText)}
      />
    </div>
  );
}

export default App;
